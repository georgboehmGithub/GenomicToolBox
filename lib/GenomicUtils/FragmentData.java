package GenomicUtils;

import GenomicAnnotations.Gene;
import GenomicAnnotations.Region;
import GenomicAnnotations.RegionVector;

public class FragmentData {
    private Gene curGene;
    private String strand;
    private RegionVector tFwRegvec;
    private RegionVector tRvRegvec;
    private Region fLocalRegion;
    private String fwRead;
    private String rvRead;
    private String fwMut;
    private String rvMut;
    private boolean invalidSequence;
    private final int readLength;

    public FragmentData(int readLength) {
        this.readLength = readLength;
        this.invalidSequence = false;
    }

    public boolean isInvalidSequence() {
        return invalidSequence;
    }

    public Gene getCurGene() {
        return curGene;
    }

    public String getRvRead() {
        return rvRead;
    }

    public String getStrand() {
        return strand;
    }

    public String getFwRead() {
        return fwRead;
    }

    public Region getfLocalRegion() {
        return fLocalRegion;
    }

    public int getReadLength() {
        return readLength;
    }

    public String getRvMut() {
        return rvMut;
    }

    public String getFwMut() {
        return fwMut;
    }

    public RegionVector gettFwRegvec() {
        return tFwRegvec;
    }

    public RegionVector gettRvRegvec() {
        return tRvRegvec;
    }

    public void setReads(String fwRead, String fwMut, String rvRead, String rvMut) {
        this.fwRead = fwRead;
        this.rvRead = rvRead;
        this.fwMut = fwMut;
        this.rvMut = rvMut;
    }

    public void setInvalidSequence(boolean invalidSequence) {
        this.invalidSequence = invalidSequence;
    }

    public void settFwRegvec(RegionVector tFwRegvec) {
        this.tFwRegvec = tFwRegvec;
    }

    public void settRvRegvec(RegionVector tRvRegvec) {
        this.tRvRegvec = tRvRegvec;
    }

    public void setfLocalRegion(Region fLocalRegion) {
        this.fLocalRegion = fLocalRegion;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public void setCurGene(Gene curGene) {
        this.curGene = curGene;
    }
}
