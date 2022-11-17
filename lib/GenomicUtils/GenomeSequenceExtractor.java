package GenomicUtils;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;

public class GenomeSequenceExtractor {
    private RandomAccessFile raf;
    IndexFile fidx;

    public GenomeSequenceExtractor(File fasta, IndexFile fidx) throws FileNotFoundException {
        raf = new RandomAccessFile(fasta, "r");
        this.fidx = fidx;
    }

    /**
     * Return fasta index entry for given chromosome
     *
     * @param chr chromosome of current gene
     * @return chromsome entry of fasta index file
     */
    public String[] getChr(String chr) {
        return fidx.getFidxEntries().get(chr);
    }

    /**
     * gtf coordinates 1-based, end inclusive. java 0-based, end exclusive.
     * start and end are supplied gtf coordinates.
     *
     * @param chr   which chr entry to look for in indexfile
     * @param start start of searched sequence (gene/transcript/exon).
     * @param end   end of searched sequence (gene/transcript/exon)
     * @return sequence
     */
    public String getSequence(String chr, long start, long end) {
        String sequence = "";
        try {
            start = start - 1; // 1 based gtf indexing

            String[] chrEntry = getChr(chr);
            long chrStart = Long.parseLong(chrEntry[2]);

            long newLinesStart = (start) / 60;
            long newLinesEnd = (end) / 60;

            long startForSeek = chrStart + start + newLinesStart;

            long byteASize = end - start + (newLinesEnd - newLinesStart);
            raf.seek(startForSeek);
            byte[] bytes = new byte[(int)byteASize];
            raf.read(bytes);
            String stringValue = new String(bytes);
            sequence = stringValue.replace("\n", "");

        } catch (IOException e) {
            e.printStackTrace();
        }

        return sequence;
    }
}
