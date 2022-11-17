package GenomicUtils;

import java.io.*;
import java.util.HashMap;
import java.util.zip.GZIPInputStream;

public class FastaParser {
    final String filePath;

    public FastaParser(String filePath) {
        this.filePath = filePath;
    }

    /**
     * Parses fasta header as key and fasta sequence as value,
     * given a fasta file into a hashmap.
     * @return Hashmap containing all fasta entries
     * @throws IOException
     */
    public HashMap<String, String> fastaToHm() throws IOException {
        var fastaHm = new HashMap<String, String>();
        BufferedReader in = new BufferedReader(new InputStreamReader(
                new GZIPInputStream(new FileInputStream(filePath))));
        String content;
        String id = "";
        String seq = "";
        while ((content = in.readLine()) != null) {
            if (content.startsWith(">")) {
                if (seq.length() > 0 && id.length() > 0) {
                    fastaHm.put(id,seq);
                }
                seq = "";
                id = content.substring(1,16);
            } else {
                seq += content;
            }

        }
        fastaHm.put(id,seq);
        return fastaHm;
    }
}
