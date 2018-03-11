import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Map;

public class CreateDictionary {

	static HashMap<Integer, Integer> segment_length = new HashMap<>();
	static int K = 200;
	static HashMap<String, Integer> output = new HashMap<String, Integer>();

	public static void main(String[] args) {
		try {
			BufferedReader in = Files.newBufferedReader(Paths.get(args[0]), StandardCharsets.UTF_8);
			BufferedWriter out = Files.newBufferedWriter(Paths.get(args[1]), StandardCharsets.UTF_8);
			K = Integer.parseInt(args[2]);
			String input = null;
			int tag = -1, len = 0;
			boolean read_segments = true;
			while ((input = in.readLine()) != null) {
				if (input.startsWith("S ")) {
					if (tag != -1) {
						segment_length.put(tag, len);
					}
					String str[] = input.split(" ");
					tag = Integer.parseInt(str[1]);
					len = str[2].length();
				} else if (read_segments) {
					len += input.length();
				}
				if (input.startsWith("P ")) {
					if (tag != -1) {
						segment_length.put(tag, len);
						tag = -1;
					}
					read_segments = false;
					String temp[] = input.split(" ");
					String path[] = temp[2].split(",");
					int path_length = path.length;

					for (int i = 0; i < path_length; i++) {
						path[i] = path[i].trim().substring(0, path[i].length() - 1);
					}

					for (int i = 0; i < path_length; i++) {
						int val = segment_length.get(Integer.parseInt(path[i]));
						StringBuilder current = new StringBuilder(path[i]);
						if (val >= K) {
							output.put(current.toString(), output.getOrDefault(current.toString(), 0) + val - K + 1);
						}
						int j = i + 1;
						while (j < path.length) {
							int sum = 0;
							for (int k = i + 1; k <= j - 1; k++) {
								sum += segment_length.get(Integer.parseInt(path[k]));
							}
							int check = sum + segment_length.get(Integer.parseInt(path[i]))
									+ segment_length.get(Integer.parseInt(path[j]));
							if (sum + 2 <= K && check >= K) {
								current.append(":" + path[j]);
								int value = Math.min(segment_length.get(Integer.parseInt(path[i])),
										segment_length.get(Integer.parseInt(path[j])));
								value = Math.min(value, check - K + 1);
								output.put(current.toString(), output.getOrDefault(current.toString(), 0) + value);
								j++;
							} else {
								current.append(":" + path[j]);
								j++;
							}
							if (sum + 2 > K) {
								break;
							}
						}
					}
				}
			}
			System.out.println();
			for (Map.Entry<String, Integer> entry : output.entrySet()) {
				out.write(entry.getKey() + "-" + entry.getValue() + "\n");
			}
			out.flush();
		} catch (Exception e) {
			e.printStackTrace(System.out);
		}
	}
}
