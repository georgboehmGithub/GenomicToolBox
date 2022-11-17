package GenomicAnnotations;

public class Transcript {

    private String id;
    private String name;
    private int start;
    private int end;
    private String protId;
    private RegionVector introns;
    private RegionVector exons;
    private CDS codingSequence;

    public Transcript(String id, String name, int start, int end) {
        this.id = id;
        this.name = name;
        this.start = start;
        this.end = end;
        this.introns = new RegionVector();
        this.exons = new RegionVector();
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    /**
     * Calculates "inverse" RegionVector object to CDS RegionVector.
     *
     * @return RegionVector of the transcript's introns
     */
    public RegionVector getIntrons() {
        if (introns.getRegions().size() == 0) {
            if (codingSequence != null) {
                codingSequence.getCdsRegions().sortRegionVector();
                Region prev = codingSequence.getCdsRegions().getRegions().get(0);
                for (int i = 1; i < codingSequence.getCdsRegions().regions.size(); i++) {
                    Region cur = codingSequence.getCdsRegions().getRegions().get(i);
                    Region newIntron = new Region(prev.getEnd() + 1, cur.getStart());
                    introns.regions.add(newIntron);
                    prev = cur;
                }
            }
        }
        return introns;
    }

    public void addExon(Region exon) {
        this.exons.getRegions().add(exon);
    }

    // Getters
    public String getProtId() {
        return protId;
    }

    public String getName() {
        return name;
    }

    public String getId() {
        return id;
    }

    public RegionVector getExons() {
        return exons;
    }

    public CDS getCDS() {
        return codingSequence;
    }

    // Setters

    public void setProtId(String protId) {
        this.protId = protId;
    }

    public void setCodingSequence(CDS codingSequence) {
        this.codingSequence = codingSequence;
    }
}
