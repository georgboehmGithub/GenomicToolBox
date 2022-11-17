package GenomicUtils;
import java.io.*;
import java.util.HashMap;

public class IndexFile {
    private HashMap<String, String[]> fidxEntries;  //key: Chr number, Value entire entry in idx file

    /**
     * Parses a Indexfile for further investigations on a genome
     * with a additional random access file.
     * @param fidxPath path of the fasta index file
     */
    public IndexFile(String fidxPath) {
        fidxEntries = new HashMap<>();
        File fidx = new File(fidxPath);
        try (BufferedReader br = new BufferedReader(new FileReader(fidx))) {
            String line;

            while ((line = br.readLine()) != null) {
                String[] entry = line.split("\t");
                fidxEntries.put(entry[0],entry);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public HashMap<String, String[]> getFidxEntries() {
        return fidxEntries;
    }
}
