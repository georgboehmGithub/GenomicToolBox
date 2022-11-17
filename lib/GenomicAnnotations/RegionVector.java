package GenomicAnnotations;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Objects;

public class RegionVector {
    ArrayList<Region> regions;
    private String strand;

    public RegionVector() {
        this.regions = new ArrayList<>();
    }

    public RegionVector(ArrayList<Region> regs) {
        this.regions = regs;
    }

    public ArrayList<Region> getRegions() {
        return regions;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    /**
     * Sorts regions in ascending order with Mergesort.
     */
    public void sortRegionVector() {
        Collections.sort(this.getRegions());
    }

    public void clear() {
        this.regions = new ArrayList<>();
    }

    public RegionVector mergedRv() {
        this.sortRegionVector();

        RegionVector copy = new RegionVector();
        for (Region re : this.getRegions()) {
            copy.getRegions().add(new Region(re.getStart(), re.getEnd()));
        }
        RegionVector result = new RegionVector();
        Region newRegion = copy.getRegions().get(0);
        result.getRegions().add(newRegion);
        for (Region re : copy.getRegions()) {
            if (newRegion.getEnd() + 1 >= re.getStart()) {
                newRegion.setEnd(Math.max(newRegion.getEnd(), re.getEnd()));
            } else {
                newRegion = re;
                result.getRegions().add(newRegion);
            }
        }
        return result;
    }

    public boolean overlapsWith(RegionVector rv) {
        for (Region re : this.regions) {
            for (Region rvRe : rv.getRegions()) {
                if (rvRe.getStart() > re.getEnd()) break;
                else if (rvRe.getStart() >= re.getStart()) return true;
                else if (rvRe.getEnd() >= re.getStart() && rvRe.getEnd() <= re.getEnd()) return true;
                else if (re.getEnd() <= rvRe.getEnd()) return true;
            }
        }
        return false;
    }

    public RegionVector cut(int start, int end) {
        Region cutRegion = new Region(start, end);
        RegionVector rvCut = new RegionVector();
        this.sortRegionVector();
        Region curReg = new Region(-1, -1);
        Region lastReg = curReg;
        int startIndex = -1;
        int cutStart = -1;
        // declare starting region
        for (int i = 0; i < this.getRegions().size(); i++) {
            lastReg = curReg;
            curReg = regions.get(i);
            if (cutRegion.overlapsWith(curReg)) {
                startIndex = i;
                cutStart = Math.max(curReg.getStart(), start);
                break;
            }
            if (lastReg.getStart() != -1 && lastReg.getEnd() < start && curReg.getStart() > end) {
                return new RegionVector();
            }
        }
        // add first region
        if (curReg.getEnd() <= end) {
            rvCut.regions.add(new Region(cutStart, curReg.getEnd()));
        } else {
            rvCut.regions.add(new Region(cutStart, end));
            return rvCut;
        }
        // add all regions until end is reached
        if (startIndex != -1) {
            for (int j = startIndex + 1; j < this.getRegions().size(); j++) {
                curReg = regions.get(j);
                if (curReg.getEnd() <= end) {
                    rvCut.regions.add(curReg);
                } else if (curReg.getEnd() >= end && curReg.getStart() <= end) {
                    rvCut.regions.add(new Region(curReg.getStart(), end));
                    break;
                } else {    // both higher
                    break;
                }
            }
        }
        return rvCut;
    }

    @Override
    public String toString() {
        String result = "";
        for (Region region : regions) {
            result += "" + region.getStart() + "-" + region.getEnd() + "|";
        }
        result = result.substring(0, result.length() - 1);
        return result;
    }

    public RegionVector getInverseRv() {
        RegionVector inverseRv = new RegionVector();
        this.sortRegionVector();
        Region prev = regions.get(0);
        for (int i = 1; i < regions.size(); i++) {
            Region cur = regions.get(i);
            Region newIntron = new Region(prev.getEnd() + 1, cur.getStart() - 1);
            inverseRv.regions.add(newIntron);
            prev = cur;
        }
        return inverseRv;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        RegionVector that = (RegionVector) o;
        return Objects.equals(regions, that.regions);
    }

    public boolean equalsRv(RegionVector rv) {
        if (this.regions.size() != rv.regions.size()) {
            return false;
        } else {
            for (int i = 0; i < regions.size(); i++) {
                if (!(regions.get(i).equals(rv.regions.get(i)))) {
                    return false;
                }
            }
        }
        return true;
    }

    @Override
    public int hashCode() {
        return Objects.hash(regions, strand);
    }

    public static void main(String[] args) {
        ArrayList<Region> a = new ArrayList<>();
        RegionVector v = new RegionVector();
        v.getRegions().add(new Region(1, 3));
        v.getRegions().add(new Region(6, 9));
        v.getRegions().add(new Region(8, 10));
        v.getRegions().add(new Region(25, 30));
        v.getRegions().add(new Region(21, 31));

        //v.mergedRv();
        RegionVector cut = v.cut(4, 24);
    }
}
