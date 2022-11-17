package Events;

import GenomicAnnotations.Gene;
import GenomicAnnotations.Region;
import GenomicAnnotations.Transcript;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.TreeSet;

public class ExonSkippingEvent {
    private final HashSet<String> SVSet;
    private final HashSet<String> WTSet;
    private final Region SVIntron;
    private final Gene gene;
    private ArrayList<Region> minSkippedExons;
    private ArrayList<Region> maxSkippedExons;
    private int minSkippedBases;
    private int maxSkippedBases;

    public ExonSkippingEvent(Gene gene, HashSet<String> SV, HashSet<String> WT, Region SVIntron) {
        this.SVSet = SV;
        this.WTSet = WT;
        this.SVIntron = SVIntron;
        this.gene = gene;
        this.minSkippedExons = new ArrayList<>();
        this.maxSkippedExons = new ArrayList<>();
        this.minSkippedBases = Integer.MAX_VALUE;
        this.maxSkippedBases = Integer.MIN_VALUE;
    }

    @Override
    public String toString() {
        return getSVIntron().toString() + "\t" + getWtIntrons() + "\t" + getWtProts() + "\t"
                + getSVProts() + "\t" + minSkippedExons.size() + "\t" + maxSkippedExons.size() + "\t"
                + getMinSkippedBases() + "\t" + getMaxSkippedBases();
    }

    public void setMinSkippedExons(ArrayList<Region> minSkippedExons) {
        this.minSkippedExons = minSkippedExons;
    }

    public void setMaxSkippedExons(ArrayList<Region> maxSkippedExons) {
        this.maxSkippedExons = maxSkippedExons;
    }

    public Gene getGene() {
        return gene;
    }

    public Region getSVIntron() {
        return SVIntron;
    }

    public ArrayList<Region> getMinSkippedExons() {
        return minSkippedExons;
    }

    public ArrayList<Region> getMaxSkippedExons() {
        return maxSkippedExons;
    }

    //TODO: Report wtsStrings, since Region Treeset seems to cause problems

    /**
     * puts together the wt introns within the SV intron, separated by '|' as
     * start:end
     *
     * @return wt introns in above format
     */
    public StringBuilder getWtIntrons() {
        var wtsStrings = new TreeSet<String>();
        for (String transcriptId : WTSet) {
            Transcript actTranscript = gene.getTranscriptById(transcriptId);
            for (Region intron : actTranscript.getIntrons().getRegions()) {
                if (intron.getStart() >= SVIntron.getStart() && intron.getEnd() <= SVIntron.getEnd()) {
                    wtsStrings.add(intron.getStart() + ":" + intron.getEnd() + "|");
                }
            }
        }
        StringBuilder wtsSb = new StringBuilder();
        for (String intron : wtsStrings) {
            wtsSb.append(intron);
        }
        wtsSb.setLength(wtsSb.length() - 1);
        return wtsSb;
    }

    /**
     * generates ids of the wt cds-s, separated by '|'
     *
     * @return Stringbuilder with wt prot ids
     */
    public StringBuilder getWtProts() {
        StringBuilder wtProts = new StringBuilder();
        for (String wt : WTSet) {
            Transcript curTrans = gene.getTranscriptById(wt);
            wtProts.append(curTrans.getProtId() + "|");
        }
        wtProts.setLength(wtProts.length() - 1);
        return wtProts;
    }

    /**
     * generates ids of the sv cds-s, separated by '|'
     *
     * @return Stringbuilder with sv prot ids
     */
    public StringBuilder getSVProts() {
        StringBuilder svProts = new StringBuilder();
        for (String sv : SVSet) {
            Transcript curTrans = gene.getTranscriptById(sv);
            svProts.append(curTrans.getProtId() + "|");
        }
        svProts.setLength(svProts.length() - 1);
        return svProts;
    }

    public void setMinSkippedBases(int minSkippedBases) {
        this.minSkippedBases = minSkippedBases;
    }

    public void setMaxSkippedBases(int maxSkippedBases) {
        this.maxSkippedBases = maxSkippedBases;
    }

    //TODO: add both to one function
    public int calcMinSkippedBases(ArrayList<Region> cand) {
        int potential = 0;
        for (Region exon : cand) {
            potential += exon.getEnd() - exon.getStart() + 1;
        }
        return potential;
    }

    public int calcMaxSkippedBases(ArrayList<Region> cand) {
        int potential = 0;
        for (Region exon : cand) {
            potential += exon.getEnd() - exon.getStart() + 1;
        }
        return potential;
    }

    public void minMaxSkippedExons(HashSet<String> WTSet, Region SVIntron) {
        ArrayList<Region> candMinExons = new ArrayList<>();
        ArrayList<Region> candMaxExons = new ArrayList<>();

        for (String id : WTSet) {
            Transcript actTranscript = gene.getTranscriptById(id);
            for (Region exon : actTranscript.getExons().getRegions()) {
                if (exon.getStart() >= SVIntron.getStart() && exon.getEnd() <= SVIntron.getEnd()) {
                    candMinExons.add(exon);
                    candMaxExons.add(exon);
                }
            }
            // Minimum
            if (getMinSkippedExons().size() == 0) {
                setMinSkippedExons(candMinExons);
            }
            if (candMinExons.size() <= getMinSkippedExons().size()) {
                setMinSkippedExons(candMinExons);
            }

            int potMinSkippedBases = calcMinSkippedBases(candMinExons);
            if (getMinSkippedBases() > potMinSkippedBases) setMinSkippedBases(potMinSkippedBases);

            // Maximum
            if (getMaxSkippedExons().size() == 0) {
                setMaxSkippedExons(candMaxExons);
            }
            if (candMaxExons.size() >= getMaxSkippedExons().size()) {
                setMaxSkippedExons(candMaxExons);
            }
            int potMaxSkippedBases = calcMaxSkippedBases(candMaxExons);
            if (getMaxSkippedBases() < potMaxSkippedBases) setMaxSkippedBases(potMaxSkippedBases);

            candMinExons = new ArrayList<Region>();
            candMaxExons = new ArrayList<Region>();
        }
    }

    public int getMinSkippedBases() {
        return minSkippedBases;
    }

    public int getMaxSkippedBases() {
        return maxSkippedBases;
    }
}
