/**
 * @file    winSketch.hpp
 * @brief   routines to index the reference 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WIN_SKETCH_HPP 
#define WIN_SKETCH_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <cassert>
#include <zlib.h>  

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/ThreadPool.hpp"
#include "map/include/string_view.hpp"

//External includes
#include "common/kseq.h"
#include "common/murmur3.h"
#include "common/prettyprint.hpp"
#include "common/sparsehash/dense_hash_map"
#include "sparsepp/spp.h"

namespace skch
{
    struct Position {
        uint32_t transcript_id_;
        uint32_t pos_;

        Position() {
            transcript_id_ = std::numeric_limits<decltype(transcript_id_)>::max();
            pos_ = std::numeric_limits<decltype(pos_)>::max();
        }

        Position(uint32_t tid, uint32_t tpos, bool torien) {
            transcript_id_ = tid;
            pos_ = tpos;
            setOrientation(torien);
        }

        void setOrientation(bool orientation) {
            if (orientation) {
                pos_ |= 1 << 31;
            } else {
                pos_ &= 0x7FFFFFFF;
            }
        }

        inline uint32_t transcript_id() { return transcript_id_; }
        inline uint32_t pos() { return (pos_ & 0x7FFFFFFF); }
        inline bool orientation() { return (pos_ & 0x80000000); }

        template <class Archive> void serialize(Archive& ar) {
            ar(transcript_id_, pos_);
        }
    };

    spp::sparse_hash_map<uint64_t, std::vector<Position>> contig2pos;

    KSEQ_INIT(gzFile, gzread)
        /**
         * @class     skch::Sketch
         * @brief     sketches and indexes the reference (subject sequence)
         * @details  
         *            1.  Minimizers are computed in streaming fashion
         *                Computing minimizers is using double ended queue which gives
         *                O(reference size) complexity
         *                Algorithm described here:
         *                https://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html
         *
         *            2.  Index hashes into appropriate format to enable fast search at L1 mapping stage
         */
        class Sketch
        {
            //private members

            //algorithm parameters
            const skch::Parameters &param;

            //Ignore top % most frequent minimizers while lookups
            const float percentageThreshold = 0.001;

            //Minimizers that occur this or more times will be ignored (computed based on percentageThreshold)
            int freqThreshold = std::numeric_limits<int>::max();

            //Make the default constructor private, non-accessible
            Sketch();

            public:

            typedef std::vector< MinimizerInfo > MI_Type;
            using MIIter_t = MI_Type::const_iterator;

            //Keep sequence length, name that appear in the sequence (for printing the mappings later)
            std::vector< ContigInfo > metadata;

            /*
             * Keep the information of what sequences come from what file#
             * Example [a, b, c] implies 
             *  file 0 contains 0 .. a-1 sequences
             *  file 1 contains a .. b-1 
             *  file 2 contains b .. c-1
             */
            std::vector< int > sequencesByFileInfo;

            //Index for fast seed lookup (unordered_map)
            /*
             * [minimizer #1] -> [pos1, pos2, pos3 ...]
             * [minimizer #2] -> [pos1, pos2...]
             * ...
             */
            using MI_Map_t = google::dense_hash_map< MinimizerMapKeyType, MinimizerMapValueType >;
            MI_Map_t minimizerPosLookupIndex,gfaLookupIndex;


            // Avoiding un-necessary stream creation + replacing strings with string view
            // is a bit > than a 2x win!
            // implementation from : https://marcoarena.wordpress.com/tag/string_view/
            std::vector<stx::string_view> split(stx::string_view str,
                    char delims) {
                std::vector<stx::string_view> ret;

                stx::string_view::size_type start = 0;
                auto pos = str.find_first_of(delims, start);
                while (pos != stx::string_view::npos) {
                    if (pos != start) {
                        ret.push_back(str.substr(start, pos - start));
                    }
                    start = pos + 1;
                    pos = str.find_first_of(delims, start);
                }
                if (start < str.length()) {
                    ret.push_back(str.substr(start, str.length() - start));
                }
                return ret;
            }

            std::vector<std::pair<uint64_t, bool>>
                explode(const stx::string_view str, const char& ch) {
                    std::string next;
                    std::vector<std::pair<uint64_t, bool>> result;
                    // For each character in the string
                    for (auto it = str.begin(); it != str.end(); it++) {
                        // If we've hit the terminal character
                        if (*it == '+' or *it == '-') {
                            bool orientation = true;
                            // If we have some characters accumulated
                            // Add them to the result vector
                            if (!next.empty()) {
                                if (*it == '-') {
                                    orientation = false;
                                }
                                result.emplace_back(std::stoll(next), orientation);
                                next.clear();
                            }
                        } else if (*it != ch) {
                            // Accumulate the next character into the sequence
                            next += *it;
                        }
                    }
                    if (!next.empty())
                        result.emplace_back(std::stoll(next),
                                true); // this case shouldn't even happen
                    return result;
                }

            struct PackedContigInfo {
                size_t fileOrder;
                size_t offset;
                uint32_t length;
            };

            private:

            /**
             * Keep list of minimizers, sequence# , their position within seq , here while parsing sequence 
             * Note : position is local within each contig
             * Hashes saved here are non-unique, ordered as they appear in the reference
             */
            MI_Type minimizerIndex,gfaIndex;

            //Frequency histogram of minimizers
            //[... ,x -> y, ...] implies y number of minimizers occur x times
            std::map<int, int> minimizerFreqHistogram;
            std::map<std::string, std::string> mymap;

            public:

            /**
             * @brief   constructor
             *          also builds, indexes the minimizer table
             */
            Sketch(const skch::Parameters &p) 
                :
                    param(p) {
                        this->build();
                        this->index();
                        this->computeFreqHist();
                    }

            private:

            /**
             * @brief     build the sketch table
             * @details   compute and save minimizers from the reference sequence(s)
             *            assuming a fixed window size
             */
            void build()
            {
                //sequence counter while parsing file
                seqno_t seqCounter = 0;

                //Create the thread pool 
                ThreadPool<InputSeqContainer, MI_Type> threadPool( [this](InputSeqContainer* e) {return buildHelper(e);}, param.threads);

                spp::sparse_hash_map<uint64_t, std::vector<std::pair<uint64_t, bool>>> path;
                spp::sparse_hash_map<uint64_t, PackedContigInfo> contigid2seq;
                std::vector<std::string> refMap;
                std::vector<uint32_t> refLengths;
                size_t k;

                for(const auto &fileName : param.gfaSequences)
                {
                    FILE *file = fopen(fileName.c_str(), "r");
                    std::string str;
                    std::ifstream input(fileName.c_str());
                    std::ofstream myfile,myfile1;
                    myfile.open ("temp.fa",ios::out);
                    myfile1.open("mapping.txt",ios::out);
                    size_t contig_cnt{0};
                    size_t ref_cnt{0};
                    size_t contig_ctr{0};
                    uint64_t maxnid{0};

                    while(getline(input,str)){
                        char firstC = str[0];
                        if (firstC != 'S' and firstC != 'P'){
                            continue;
                        }
                        stx::string_view lnview(str);
                        std::vector<stx::string_view> splited = split(lnview, '\t');
                        string tag = splited[0].to_string();
                        string id = splited[1].to_string();
                        string value = splited[2].to_string();

                        /*kseq_t *seq1 = (kseq_t*)calloc(1, sizeof(kseq_t));
                          kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));    
                          gzFile fp = gzdopen(fileno(file), "r");
                          ks->f = fp;
                          ks->buf = (char*)malloc(1024);
                          seq1->f = ks;*/

                        if(tag=="S"){
                            mymap[id]=value;
                            myfile << ">" << id << "\n" << value << endl;
                            auto nid = std::stoull(splited[1].to_string());
                            (void)nid;
                            auto clen = splited[2].length();
                            contigid2seq[nid] = {contig_ctr, 0, static_cast<uint32_t>(clen)};
                            ++contig_ctr;
                        }
                        else if(tag=="P"){
                            k=param.kmerSize;
                            auto pvalue = splited[2];
                            std::vector<std::pair<uint64_t, bool>> contigVec = explode(pvalue, ',');
                            path[ref_cnt] = contigVec;
                            uint32_t refLength{0};
                            bool firstContig{true};

                            for (auto& ctig : contigVec) {
                                int32_t l = contigid2seq[ctig.first].length - (firstContig ? 0 : (k-1));
                                refLength += l;
                                firstContig = false;
                            }

                            refLengths.push_back(refLength);
                            refMap.push_back(id);
                            ref_cnt++;
                        }
                    }

                    //Make the pufferfish mapping in file.
                    uint64_t pos = 0;
                    uint64_t accumPos;
                    uint64_t currContigLength = 0;
                    uint64_t total_output_lines = 0;
                    for (auto const& ent : path) {
                        const uint64_t& tr = ent.first;
                        const std::vector<std::pair<uint64_t, bool>>& contigs = ent.second;
                        accumPos = 0;
                        for (size_t i = 0; i < contigs.size(); i++) {
                            if (contig2pos.find(contigs[i].first) == contig2pos.end()) {
                                contig2pos[contigs[i].first] = {};
                                total_output_lines += 1;
                            }
                            if (contigid2seq.find(contigs[i].first) == contigid2seq.end()) {
                                std::cerr << contigs[i].first << "\n";
                            }
                            pos = accumPos;
                            currContigLength = contigid2seq[contigs[i].first].length;
                            accumPos += currContigLength - k;
                            (contig2pos[contigs[i].first])
                                .push_back(Position(tr, pos, contigs[i].second));
                        }
                    }
                    myfile1<<"Contig ID : space_seperated(transcript_id,position)"<<endl;
                    for(auto &e : contig2pos){
                        std::vector<Position> temp=e.second;
                        myfile1 << e.first<<":";
                        for(auto &v :temp){
                            myfile1 << v.transcript_id()<<","<<v.pos()<<" ";
                        }
                        myfile1 << endl;
                    }

                    myfile.close();
                    myfile1.close();
                }

                for(const auto &fileName : param.refSequences)
                {

#ifdef DEBUG
                    std::cout << "INFO, skch::Sketch::build, building minimizer index for " << fileName << std::endl;
#endif

                    //Open the file using kseq
                    //FILE *file = fopen("temp.fa", "r");
                    //gzFile fp = gzdopen(fileno(file), "r");
                    //kseq_t *seq = kseq_init(fp);


                    //size of sequence
                    offset_t len;

                    for(auto &e:mymap) 
                    {
                        len=e.second.length();
                        //Save the sequence name
                        metadata.push_back( ContigInfo{e.first, (offset_t)len} );

                        //Is the sequence too short?
                        if(len < param.windowSize || len < param.kmerSize)
                        {
#ifdef DEBUG
                            //cout<<len<<":"<<param.windowSize<<":"<<param.kmerSize<<":"<<seq->seq.s<<endl;
                            std::cout << "WARNING, skch::Sketch::build, found an unusually short sequence relative to kmer and window size" << std::endl;
#endif
                            seqCounter++;
                            continue;  
                        }
                        else
                        {
                            threadPool.runWhenThreadAvailable(new InputSeqContainer(e.second.c_str(), e.first.c_str(), len, seqCounter));

                            //Collect output if available
                            while ( threadPool.outputAvailable() )
                                this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());
                        }

                        seqCounter++;
                    }

                    sequencesByFileInfo.push_back(seqCounter);

                    //kseq_destroy(seq);  
                    //gzclose(fp); //close the file handler 
                    //fclose(file);
                }


                //Collect remaining output objects
                while ( threadPool.running() )
                    this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());

                std::cout << "INFO, skch::Sketch::build, minimizers picked from reference = " << minimizerIndex.size() << std::endl;

            }

            /**
             * @brief               function to compute minimizers given input sequence object
             * @details             this function is run in parallel by multiple threads
             * @param[in]   input   input read details
             * @return              output object containing the mappings
             */
            MI_Type* buildHelper(InputSeqContainer *input)
            {
                MI_Type* thread_output = new MI_Type();

                //Compute minimizers in reference sequence
                skch::CommonFunc::addMinimizers(*thread_output, &(input->seq[0u]), input->len, param.kmerSize, param.windowSize, param.alphabetSize, input->seqCounter);

                return thread_output;
            }

            /**
             * @brief                 routine to handle thread's local minimizer index
             * @param[in] output      thread local minimizer output
             */
            void buildHandleThreadOutput(MI_Type* output)
            {
                this->minimizerIndex.insert(this->minimizerIndex.end(), output->begin(), output->end());
                delete output;
            }

            /**
             * @brief   build the index for fast lookups using minimizer table
             */
            void index()
            {
                //Parse all the minimizers and push into the map
                minimizerPosLookupIndex.set_empty_key(0);

                for(auto &e : minimizerIndex)
                {
                    // [hash value -> info about minimizer]
                    minimizerPosLookupIndex[e.hash].push_back( 
                            MinimizerMetaData{e.seqId, e.wpos, e.strand});
                }

                std::cout << "INFO, skch::Sketch::index, unique minimizers = " << minimizerPosLookupIndex.size() << std::endl;

                //Writing mashmap index to a file.
                std::ofstream myfile;
                myfile.open ("mashmap_dict.txt",ios::out);

                myfile << "hash:space_seperated_list[seqId,wpos,strand]" <<endl;
                for(auto &e : this->minimizerPosLookupIndex){
                    myfile << e.first<<":";
                    for(auto &f : e.second){
                        myfile << f.seqId<<","<<f.wpos<<","<<f.strand<<" ";
                    }
                    myfile << endl;
                }
                myfile.close();
            }

            /**
             * @brief   report the frequency histogram of minimizers using position lookup index
             *          and compute which high frequency minimizers to ignore
             */
            void computeFreqHist()
            {

                //1. Compute histogram

                for(auto &e : this->minimizerPosLookupIndex)
                    this->minimizerFreqHistogram[e.second.size()] += 1;

                std::cout << "INFO, skch::Sketch::computeFreqHist, Frequency histogram of minimizers = " <<  *this->minimizerFreqHistogram.begin() <<  " ... " << *this->minimizerFreqHistogram.rbegin() << std::endl;

                //2. Compute frequency threshold to ignore most frequent minimizers

                int64_t totalUniqueMinimizers = this->minimizerPosLookupIndex.size();
                int64_t minimizerToIgnore = totalUniqueMinimizers * percentageThreshold / 100;

                int64_t sum = 0;

                //Iterate from highest frequent minimizers
                for(auto it = this->minimizerFreqHistogram.rbegin(); it != this->minimizerFreqHistogram.rend(); it++)
                {
                    sum += it->second; //add frequency
                    if(sum < minimizerToIgnore)
                    {
                        this->freqThreshold = it->first;
                        //continue
                    }
                    else if(sum == minimizerToIgnore)
                    {
                        this->freqThreshold = it->first;
                        break;
                    }
                    else
                    {
                        break;
                    }
                }

                if(this->freqThreshold != std::numeric_limits<int>::max())
                    std::cout << "INFO, skch::Sketch::computeFreqHist, With threshold " << this->percentageThreshold << "\%, ignore minimizers occurring >= " << this->freqThreshold << " times during lookup." << std::endl;
                else
                    std::cout << "INFO, skch::Sketch::computeFreqHist, With threshold " << this->percentageThreshold << "\%, consider all minimizers during lookup." << std::endl;

            }

            public:

            /**
             * @brief               search hash associated with given position inside the index
             * @details             if MIIter_t iter is returned, than *iter's wpos >= winpos
             * @param[in]   seqId
             * @param[in]   winpos
             * @return              iterator to the minimizer in the index
             */
            MIIter_t searchIndex(seqno_t seqId, offset_t winpos) const
            {
                std::pair<seqno_t, offset_t> searchPosInfo(seqId, winpos);

                /*
                 * std::lower_bound --  Returns an iterator pointing to the first element in the range
                 *                      that is not less than (i.e. greater or equal to) value.
                 */
                MIIter_t iter = std::lower_bound(this->minimizerIndex.begin(), this->minimizerIndex.end(), searchPosInfo, cmp);

                return iter;
            }

            /**
             * @brief                 check if iterator points to index end
             * @param[in]   iterator
             * @return                boolean value
             */
            bool isMinimizerIndexEnd(const MIIter_t &it) const
            {
                return it == this->minimizerIndex.end();
            }

            /**
             * @brief     Return end iterator on minimizerIndex
             */
            MIIter_t getMinimizerIndexEnd() const
            {
                return this->minimizerIndex.end();
            }

            int getFreqThreshold() const
            {
                return this->freqThreshold;
            }

            private:

            /**
             * @brief     functor for comparing minimizers by their position in minimizerIndex
             * @details   used for locating minimizers with the required positional information
             */
            struct compareMinimizersByPos
            {
                typedef std::pair<seqno_t, offset_t> P;

                bool operator() (const MinimizerInfo &m, const P &val)
                {
                    return ( P(m.seqId, m.wpos) < val);
                }

                bool operator() (const P &val, const MinimizerInfo &m)
                {
                    return (val < P(m.seqId, m.wpos) );
                }
            } cmp;

        }; //End of class Sketch
} //End of namespace skch

#endif
