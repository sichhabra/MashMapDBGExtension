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

KSEQ_INIT(gzFile, gzread)

    namespace skch
{
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

            for(const auto &fileName : param.gfaSequences)
            {
                FILE *file = fopen(fileName.c_str(), "r");
                std::string str;
                std::ifstream input(fileName.c_str());
                std::ofstream myfile;
                myfile.open ("temp.fa",ios::out);
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
                    
                    if(tag=="S"){
                        myfile << ">" << id << "\n" << value << endl;
                    }
                }
                myfile.close();
            }

            for(const auto &fileName : param.refSequences)
            {

#ifdef DEBUG
                std::cout << "INFO, skch::Sketch::build, building minimizer index for " << fileName << std::endl;
#endif

                //Open the file using kseq
                FILE *file = fopen(fileName.c_str(), "r");
                gzFile fp = gzdopen(fileno(file), "r");
                kseq_t *seq = kseq_init(fp);


                //size of sequence
                offset_t len;

                while ((len = kseq_read(seq)) >= 0) 
                {
                    //Save the sequence name
                    metadata.push_back( ContigInfo{seq->name.s, (offset_t)seq->seq.l} );

                    //Is the sequence too short?
                    if(len < param.windowSize || len < param.kmerSize)
                    {
#ifdef DEBUG
                        cout<<len<<":"<<param.windowSize<<":"<<param.kmerSize<<":"<<seq->seq.s<<endl;
                        std::cout << "WARNING, skch::Sketch::build, found an unusually short sequence relative to kmer and window size" << std::endl;
#endif
                        seqCounter++;
                        continue;  
                    }
                    else
                    {
                        threadPool.runWhenThreadAvailable(new InputSeqContainer(seq->seq.s, seq->name.s, len, seqCounter));

                        //Collect output if available
                        while ( threadPool.outputAvailable() )
                            this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());
                    }

                    seqCounter++;
                }

                sequencesByFileInfo.push_back(seqCounter);

                kseq_destroy(seq);  
                gzclose(fp); //close the file handler 
                fclose(file);
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

            /*for(auto &e : this->minimizerPosLookupIndex){
              cout<< e.first<<":";
              for(auto &f : e.second){
              cout<< f.seqId<<","<<f.wpos<<","<<f.strand<<" ";
              }
              cout<<endl;
              }
             */
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
