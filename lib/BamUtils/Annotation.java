package BamUtils;

import GenomicAnnotations.Gene;
import GenomicAnnotations.Region;
import GenomicAnnotations.RegionVector;
import GenomicAnnotations.Transcript;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;
import java.util.ArrayList;
import java.util.HashMap;

public class Annotation {
    private SAMRecord firstOfPair;
    private SAMRecord secondOfPair;
    private RegionVector fopRv;
    private RegionVector sopRv;
    private String strandness;
    private boolean antisenseMatch;
    private boolean intergenic;
    private HashMap<Integer, ArrayList<Gene>> relGenes;
    private HashMap<String, ArrayList<Transcript>> relTranscripts;
    private int highestFeature;
    private int cgenes;
    private int igenes;
    private int clipped;
    private int mismatches;
    private int gDist;
    private int splitCount;

    public Annotation() {
        this.strandness = "";
        this.intergenic = false;
        this.fopRv = new RegionVector();
        this.sopRv = new RegionVector();
        this.clipped = 0;
        this.mismatches = 0;
        this.splitCount = 0;
        this.cgenes = 0;
        this.igenes = 0;
        this.antisenseMatch = false;
        this.firstOfPair = null;
        this.secondOfPair = null;
        this.highestFeature = 0;
    }

    public Annotation clear() {
        return new Annotation();
    }

    /**
     * Runtime O(n) if read region vector not yet initialized
     *
     * @param fop   specifies wether first or second read of pair
     * @return      regions of alignment for read
     */
    public RegionVector getReadRv(boolean fop) {
        RegionVector readRv = fop ? fopRv : sopRv;
        SAMRecord read = fop ? firstOfPair : secondOfPair;

        if (readRv.getRegions().size() > 0) {    // rv initialized
            return readRv;
        } else {
            RegionVector readRegions = new RegionVector();
            for (AlignmentBlock ab : read.getAlignmentBlocks()) {
                readRegions.getRegions().add(new Region(ab.getReferenceStart(), ab.getReferenceStart() + ab.getLength() - 1));
            }
            readRegions = readRegions.mergedRv();
            return readRegions;
        }
    }

    /**
     * Produces annotation for split inconsistent read pair
     *
     * @return  string representing split inconsistent read pair
     */
    public String toStringSplitInconsistent() {
        return firstOfPair.getReadName() + "\t" + "split-inconsistent:true" + "\n";
    }

    /**
     * Produces annotation for intergenic read pair
     *
     * @return  string representing split inconsistent read pair
     */
    public String toStringIntergenic(int pcrIndex) {
        return firstOfPair.getReadName() + "\t" + "mm:" +
                mismatches + "\t" + "clipping:" + clipped +
                "\t" + "nsplit:" + splitCount + "\t" +
                "gcount:" + cgenes + "\t" + "gdist:" +
                gDist + "\t" + "antisense:" + antisenseMatch + "\t" + "pcrindex:" + pcrIndex + "\n";
    }

    public String toStringGenic() {
        return null;
    }

    public HashMap<String, ArrayList<Transcript>> getRelTranscripts() {
        return relTranscripts;
    }

    public void setRelTranscripts(HashMap<String, ArrayList<Transcript>> relTranscripts) {
        this.relTranscripts = relTranscripts;
    }

    public HashMap<Integer, ArrayList<Gene>> getRelGenes() {
        return relGenes;
    }

    public void setRelGenes(HashMap<Integer, ArrayList<Gene>> relGenes) {
        this.relGenes = relGenes;
    }

    public int getCgenes() {
        return cgenes;
    }

    public int getHighestFeature() {
        return highestFeature;
    }

    public void setHighestFeature(int highestFeature) {
        this.highestFeature = highestFeature;
    }

    public boolean isIntergenic() {
        return intergenic;
    }

    public void setIntergenic(boolean intergenic) {
        this.intergenic = intergenic;
    }

    public boolean isAntisenseMatch() {
        return antisenseMatch;
    }

    public void setAntisenseMatch(boolean antisenseMatch) {
        this.antisenseMatch = antisenseMatch;
    }

    public void incCGenes(int cGenesCount) {
        this.cgenes += cGenesCount;
    }

    public void incIGenes(int iGenesCount) {
        this.igenes += iGenesCount;
    }

    public int getClipped() {
        return clipped;
    }

    public void setClipped(int clipped) {
        this.clipped = clipped;
    }

    public int getMismatches() {
        return mismatches;
    }

    public void setMismatches(int mismatches) {
        this.mismatches = mismatches;
    }

    public int getgDist() {
        return gDist;
    }

    public void setgDist(int gDist) {
        this.gDist = gDist;
    }

    public int getSplitCount() {
        return splitCount;
    }

    public void setSplitCount(int splitCount) {
        this.splitCount = splitCount;
    }

    public SAMRecord getFirstOfPair() {
        return firstOfPair;
    }

    public void setFirstOfPair(SAMRecord firstOfPair) {
        this.firstOfPair = firstOfPair;
    }

    public SAMRecord getSecondOfPair() {
        return secondOfPair;
    }

    public void setSecondOfPair(SAMRecord secondOfPair) {
        this.secondOfPair = secondOfPair;
    }

    public void setStrandness(String strandness) {
        this.strandness = strandness;
    }
}
