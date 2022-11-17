package GenomicAnnotations;

import java.util.HashMap;
import java.util.HashSet;

public class Gene {
    private final String chromosome;
    private final String id;
    private final String name;
    private final String biotype;
    private final String strand;
    private int nprots;
    private final int start;
    private final int end;
    private HashMap<String, Transcript> transcripts;
    private HashMap<String, CDS> CDSs;

    public Gene(String chromosome, String id, String name, String biotype, String strand, int start, int end) {
        this.chromosome = chromosome;
        this.id = id;
        this.name = name;
        this.biotype = biotype;
        this.strand = strand;
        this.transcripts = new HashMap<>();
        this.CDSs = new HashMap<>();
        this.nprots = 0;
        this.start = start;
        this.end = end;
    }

    /**
     * Captures unique introns of a genes' set of transcripts
     *
     * @return Hashset of unique introns
     */
    public HashSet<Region> uniqueIntrons() {
        HashSet<Region> uniqueIntrons = new HashSet<>();
        HashMap<String, Transcript> geneTranscripts = transcripts;
        geneTranscripts.forEach((transcriptId, curTrans) -> {
            Transcript t = curTrans;
            t.getIntrons();
            uniqueIntrons.addAll(t.getIntrons().getRegions());
        });
        return uniqueIntrons;
    }

    @Override
    public String toString() {
        return getId() + "\t" + getName() + "\t" + getChromosome() + "\t" + getStrand() + "\t"
                + getNprots() + "\t" + getTranscripts().size() + "\t";
    }

    public String getBiotype() {
        return biotype;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public Transcript getTranscriptById(String id) {
        return transcripts.get(id);
    }

    public CDS getCdsById(String id) {
        return CDSs.get(id);
    }

    public HashMap<String, CDS> getCDSs() {
        return CDSs;
    }

    public int getNprots() {
        return nprots;
    }

    public void setNprots(int nprots) {
        this.nprots = nprots;
    }

    public String getStrand() {
        return strand;
    }

    public String getChromosome() {
        return chromosome;
    }

    public String getName() {
        return name;
    }

    public String getId() {
        return id;
    }

    public HashMap<String, Transcript> getTranscripts() {
        return transcripts;
    }
}
