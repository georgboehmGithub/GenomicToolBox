package GenomicAnnotations;

public class CDS {

    private final String proteinId;
    private int start;
    private int end;
    private final String transcriptId;
    private RegionVector cdsRegions;

    public CDS(String proteinId, String transcriptId) {
        this.proteinId = proteinId;
        this.transcriptId = transcriptId;
        this.cdsRegions = new RegionVector();
    }

    public RegionVector getCdsRegions() {
        return cdsRegions;
    }

    public void addCds(Region cds) {
        cdsRegions.getRegions().add(cds);
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }
}
