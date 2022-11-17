package GenomicAnnotations;

import augmentedTree.Interval;

import java.util.Objects;

public class Region implements Comparable<Region>, Interval {

    private int start;
    private int end;
    private int length;
    private String geneId;
    private String strand;

    public Region(int start, int end) {
        this.start = start;
        this.end = end;
        this.length = getLength();
        checkIntegrity();
    }

    public Region(int start, int end, String geneId) {
        this.start = start;
        this.end = end;
        this.length = getLength();
        this.geneId = geneId;
        checkIntegrity();
    }

    /**
     * Makes sure that regions contain logic start and end coordinates
     */
    public void checkIntegrity() {
        if (start <= end) {
            return;
        }
        System.out.println("Invalid Region, end should be higher than start: " + start + " " + end);
    }

    public String getGeneId() {
        return geneId;
    }

    public void setGeneId(String geneId) {
        this.geneId = geneId;
    }

    public boolean overlapsWith(Region re) {
        if (this.getStart() > re.getEnd()) return false;
        else if (this.getStart() >= re.getStart()) return true;
        else if (this.getEnd() >= re.getStart() && this.getEnd() <= re.getEnd()) return true;
        else if (re.getEnd() <= this.getEnd()) return true;
        return false;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int getStop() {
        return end;
    }

    public int getLength() {
        return this.getEnd() - this.getStart() + 1;
    }

    public String getStrand() {
        return strand;
    }

    public void setStart(int start) {
        this.start = start;
    }

    public void setEnd(int end) {
        this.end = end;
    }

    @Override
    public String toString() {
        return start + ":" + end;
    }

    @Override
    public int compareTo(Region o) {
        return Integer.compare(this.start, o.start);
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Region region = (Region) o;
        return start == region.start &&
                end == region.end;
    }

    @Override
    public int hashCode() {
        return Objects.hash(start, end);
    }
}
