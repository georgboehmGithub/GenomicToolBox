package BamUtils;

import GenomicAnnotations.*;
import augmentedTree.IntervalTree;
import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.SAMRecord;

import java.util.*;

public class GeneralUtils {
    private Annotation annotation;
    private GenomeAnnotations gAnnotations;
    private HashMap<RegionVector, Integer> pcrCounts;

    public GeneralUtils(GenomeAnnotations gAnnotations) {
        this.annotation = new Annotation();
        this.gAnnotations = gAnnotations;
        this.pcrCounts = new HashMap<>();
    }

    public HashMap<RegionVector, Integer> getPcrCounts() {
        return pcrCounts;
    }

    // Annotation start //
    public Annotation getAnnotation() {
        return annotation;
    }

    public void setAnnotation(Annotation newAnnot) {
        this.annotation = newAnnot;
    }

    public void setAnnotation(SAMRecord firstOfPair, SAMRecord secondOfPair, String strandness) {
        annotation.setFirstOfPair(firstOfPair);
        annotation.setSecondOfPair(secondOfPair);
        annotation.setStrandness(strandness);
    }
    // Annotation end //

    // Clipping start //
    public void readPairClipping(SAMRecord sr, SAMRecord mate) {
        annotation.setClipped(clipping(sr) + clipping(mate));
    }

    public int clipping(SAMRecord sr) {
        int startDiff = Math.abs(sr.getAlignmentStart() - sr.getUnclippedStart());
        int endDiff = Math.abs(sr.getAlignmentEnd() - sr.getUnclippedEnd());
        return startDiff + endDiff;
    }
    // Clipping end //

    // Mismatch start //
    public void readpairMismatchCount(SAMRecord sr, SAMRecord mate) {
        annotation.setMismatches(readMismatchCount(sr) + readMismatchCount(mate));
    }

    public int readMismatchCount(SAMRecord sr) {
        int nm = -1;
        nm = (Integer) sr.getAttribute("NM");
        nm = (nm != -1) ? nm : (Integer) sr.getAttribute("nM");
        nm = (nm != -1) ? nm : (Integer) sr.getAttribute("XM");
        return (nm != -1) ? nm : 0;
    }
    // Mismatch end //

    // Split count start //
    public int splitCount(RegionVector srIntrons, RegionVector mateIntrons) {
        var uniqueIntrons = new HashSet<>();
        uniqueIntrons.addAll(srIntrons.getRegions());
        uniqueIntrons.addAll(mateIntrons.getRegions());
        return uniqueIntrons.size();
    }
    // Split count end //

    // Valid read pair start //

    /**
     * Checks validity for single read
     *
     * @param sr read being checked
     * @return true: valid read, up for further processing
     */
    public boolean validRead(SAMRecord sr) {
        if (sr.isSecondaryOrSupplementary()) return false;
        if (sr.getReadUnmappedFlag()) return false;
        if (sr.getMateUnmappedFlag()) return false;
        // same strand
        if (sr.getReadNegativeStrandFlag() && sr.getMateNegativeStrandFlag() ||
                !sr.getReadNegativeStrandFlag() && !sr.getMateNegativeStrandFlag()) return false;
        return true;
    }

    /**
     * Checks for two single valid reads if their corresponding read pair is valid
     *
     * @param intervalTree     sense strand intervals of gene the read lies on
     * @param intervalTreeAnti anti sense strand intervals of read
     * @param unspecStrandness true: no strandness for read pair given
     * @return true: valid read pair, up for further processing
     */
    public boolean validReadPair(IntervalTree<Region> intervalTree, IntervalTree<Region> intervalTreeAnti, boolean unspecStrandness) {
        SAMRecord firstOfPair = annotation.getFirstOfPair();
        SAMRecord secondOfPair = annotation.getSecondOfPair();

        ArrayList<Region> igenes = new ArrayList<>();
        ArrayList<Region> cgenes = new ArrayList<>();
        int start = Math.min(Math.min(firstOfPair.getStart(), secondOfPair.getStart()), Math.min(firstOfPair.getEnd(), secondOfPair.getEnd()));
        int end = Math.max(Math.max(firstOfPair.getEnd(), secondOfPair.getEnd()), Math.max(firstOfPair.getStart(), secondOfPair.getStart()));

        igenes = intervalTree.getIntervalsSpannedBy(start, end, igenes);
        cgenes = intervalTree.getIntervalsSpanning(start, end, cgenes);
        if (unspecStrandness) { // add antisene genes
            igenes.addAll(intervalTreeAnti.getIntervalsSpannedBy(start, end, igenes));
            cgenes.addAll(intervalTreeAnti.getIntervalsSpanning(start, end, cgenes));
        }
        annotation.incCGenes(cgenes.size());
        annotation.incIGenes(igenes.size());

        if (igenes.size() > 0 && cgenes.size() == 0) {  // no valid read pair
            return false;
        } else if (igenes.size() == 0 && cgenes.size() == 0) {  // intergenic read pair
            handleIntergenic(unspecStrandness, firstOfPair, secondOfPair, intervalTree, intervalTreeAnti, start, end);
        } else {    // genic read pair
            annotation.setIntergenic(false);
            calcRelGenes(cgenes, firstOfPair, secondOfPair);
        }
        // check for antisense situation
        handleAntisense(unspecStrandness, intervalTreeAnti, start, end);
        return true;
    }
    // Valid read pair end //

    public void calcRelGenes(ArrayList<Region> cgenes, SAMRecord firstOfPair, SAMRecord secondOfPair) {
        // TODO: Sub functions "check transcriptomic" etc. for better oversight in e.g. Annotation class
        var relGenes = new HashMap<Integer, ArrayList<Gene>>();
        var relTranscripts = new HashMap<String, ArrayList<Transcript>>();
        HashSet<Transcript> subsTogether = new HashSet<>();

        var geneList = new HashSet<String>();
        for (Region re : cgenes) {
            geneList.add(re.getGeneId());
        }
        Gene curGene;
        for (String geneId : geneList) {
            var transcriptList = new HashSet<Transcript>();
            curGene = gAnnotations.getGeneById(geneId);
            for (Map.Entry<String, Transcript> transcriptEntry : curGene.getTranscripts().entrySet()) {
                transcriptList.add(transcriptEntry.getValue());
            }
            RegionVector fopRv = annotation.getReadRv(true);
            RegionVector sopRv = annotation.getReadRv(false);

            for (Transcript t : transcriptList) {
                RegionVector fopTransCut = t.getExons().cut(firstOfPair.getAlignmentStart(), firstOfPair.getAlignmentEnd());
                RegionVector sopTransCut = t.getExons().cut(secondOfPair.getAlignmentStart(), secondOfPair.getAlignmentEnd());
                if (fopTransCut.equalsRv(fopRv) && sopTransCut.equalsRv(sopRv)) {
                    subsTogether.add(t);
                }
            }

            if (subsTogether.size() > 0) { // transcriptomic
                // transcriptomic read
                annotation.setHighestFeature(2);
                if (!relGenes.containsKey(2)) {
                    relGenes.put(2, new ArrayList<>());
                }
                relGenes.get(2).add(curGene);
                relTranscripts.put(geneId, new ArrayList<>());
                relTranscripts.get(geneId).addAll(subsTogether);

            } else {    // no transcriptomic read
                int containedRegions = 0;
                var mergedRegs = new ArrayList<Region>();
                mergedRegs.addAll(fopRv.getRegions());
                mergedRegs.addAll(sopRv.getRegions());
                RegionVector mergedReadRv = new RegionVector(mergedRegs);
                mergedReadRv = mergedReadRv.mergedRv();

                IntervalTree<Region> allExons = new IntervalTree<>();
                RegionVector intervalsMerged = new RegionVector();

                for (Transcript t : transcriptList) {
                    intervalsMerged.getRegions().addAll(t.getExons().getRegions());
                }
                intervalsMerged = intervalsMerged.mergedRv();
                allExons.addAll(intervalsMerged.getRegions());

                for (Region reg : mergedReadRv.getRegions()) {
                    ArrayList<Region> overlaps = new ArrayList<>();
                    overlaps = allExons.getIntervalsSpanning(reg.getStart(), reg.getEnd(), overlaps);
                    if (overlaps.size() == 1) {  // region contained
                        containedRegions++;
                    } else {
                        break;
                    }
                }
                if (containedRegions == mergedReadRv.getRegions().size()) { // merged
                    if (!relGenes.containsKey(2)) annotation.setHighestFeature(1);
                    if (!relGenes.containsKey(1)) {
                        relGenes.put(1, new ArrayList<>());
                    }
                    relGenes.get(1).add(curGene);
                } else {    // intronic
                    if (!relGenes.containsKey(0)) {
                        relGenes.put(0, new ArrayList<>());
                    }
                    relGenes.get(0).add(curGene);
                }
            }
        }
        annotation.setRelGenes(relGenes);
        annotation.setRelTranscripts(relTranscripts);
    }

    public int calculatePcrIndices(SAMRecord firstOfPair, SAMRecord secondOfPair, String strandness) {
        RegionVector merged = new RegionVector();
        for (AlignmentBlock ab : firstOfPair.getAlignmentBlocks()) {
            merged.getRegions().add(new Region(ab.getReferenceStart(), ab.getReferenceStart() + ab.getLength() - 1));
        }
        for (AlignmentBlock ab : secondOfPair.getAlignmentBlocks()) {
            merged.getRegions().add(new Region(ab.getReferenceStart(), ab.getReferenceStart() + ab.getLength() - 1));
        }
        merged = merged.mergedRv();
        merged.setStrand(strandness);
        if (!pcrCounts.containsKey(merged)) {
            pcrCounts.put(merged, 0);
        } else {
            int curCount = pcrCounts.get(merged);
            pcrCounts.put(merged, curCount + 1);
        }
        return pcrCounts.get(merged);
    }

    /**
     * Decides for a given read pair of SAMRecords if they are split consistent
     *
     * @param sr   first part of read pair
     * @param mate second part of read pair
     * @return true: read pair is split consistent. false: no split consistency
     */
    public boolean splitConsistent(SAMRecord sr, SAMRecord mate) {
        int overlapStart = Math.max(sr.getAlignmentStart(), mate.getAlignmentStart());
        int overlapEnd = Math.min(sr.getAlignmentEnd(), mate.getAlignmentEnd());
        IntervalTree<Region> fopIntervals;
        IntervalTree<Region> sopIntervals;

        RegionVector fopRegs = annotation.getReadRv(true);
        RegionVector sopRegs = annotation.getReadRv(false);
        RegionVector inverseFopRegs = fopRegs.getInverseRv();
        RegionVector inverseSopRegs = sopRegs.getInverseRv();

        // possibly split inconsistent
        if (overlapEnd > overlapStart) {
            fopIntervals = new IntervalTree<>(fopRegs.getRegions());
            sopIntervals = new IntervalTree<>(sopRegs.getRegions());
            var overlaps = new ArrayList<Region>();

            for (Region intron : inverseFopRegs.getRegions()) {
                overlaps = sopIntervals.getIntervalsIntersecting(intron.getStart(), intron.getEnd(), overlaps);
                if (overlaps.size() > 0) {
                    return false;
                }
            }
            for (Region intron : inverseSopRegs.getRegions()) {
                overlaps = fopIntervals.getIntervalsIntersecting(intron.getStart(), intron.getEnd(), overlaps);
                if (overlaps.size() > 0) {
                    return false;
                }
            }
        }
        annotation.setSplitCount(splitCount(inverseFopRegs, inverseSopRegs));
        return true;
    }

    /**
     * Calculates the minimum distance of a SAMRecord read pair to a neighbour gene
     * on the same strand
     *
     * @param firstOfPair  first part of read pair
     * @param secondOfPair second part of read pair
     * @param it           interval tree on the same strand
     * @param start        leftmost part of the read pair
     * @param end          rightmost part of the read pair
     * @return distance to closest gene on the same strand
     */
    public int calcDistancesNeighbours(SAMRecord firstOfPair, SAMRecord secondOfPair, IntervalTree<Region> it, int start, int end) {
        ArrayList<Region> rightNeighbour = new ArrayList<>();
        ArrayList<Region> leftNeighbour = new ArrayList<>();
        rightNeighbour = it.getIntervalsRightNeighbor(start, end, rightNeighbour);
        leftNeighbour = it.getIntervalsLeftNeighbor(start, end, leftNeighbour);

        int rightMostCoord = Math.max(firstOfPair.getAlignmentEnd(), secondOfPair.getAlignmentEnd());
        int leftMostCoord = Math.min(firstOfPair.getAlignmentStart(), secondOfPair.getAlignmentStart());

        int closestRight = rightNeighbour.size() > 0 ? rightNeighbour.get(0).getStart() : -1;
        int closestLeft = leftNeighbour.size() > 0 ? leftNeighbour.get(0).getEnd() : -1;

        if (closestRight != -1) {
            for (Region re : rightNeighbour) {
                closestRight = Math.min(closestRight, re.getStart());
            }
        }
        if (closestLeft != -1) {
            for (Region re : leftNeighbour) {
                closestRight = Math.max(closestLeft, re.getEnd());
            }
        }
        int leftDist = Math.abs(leftMostCoord - closestLeft + 1);
        int rightDist = Math.abs(closestRight - rightMostCoord + 1);
        return Math.min(leftDist, rightDist);
    }

    /**
     * Decides wether readpair is contained in a gene on antisense strand
     * on the same chromosome
     *
     * @param unspecStrandness true:no strandness given
     * @param itAnti           chromosome intervals on anti sense
     * @param start            alignment start
     * @param end              alignment end
     */
    public void handleAntisense(boolean unspecStrandness, IntervalTree<Region> itAnti,
                                int start, int end) {
        HashSet<Region> cgenes = new HashSet<>();
        if (!unspecStrandness) {
            cgenes = itAnti.getIntervalsSpanning(start, end, cgenes);
            if (cgenes.size() > 0) {
                annotation.setAntisenseMatch(true);
            }
        }
    }

    /**
     * @param unspecStrandness true:no strandness given
     * @param firstOfPair      first read of readpair
     * @param secondOfPair     second read of readpair
     * @param it               sense chromosome interval tree
     * @param antiIt           antisense chromosome interval tree
     * @param start            alignment start
     * @param end              alignment end
     */
    public void handleIntergenic(boolean unspecStrandness, SAMRecord firstOfPair, SAMRecord secondOfPair,
                                 IntervalTree<Region> it, IntervalTree<Region> antiIt, int start, int end) {
        // intergenic
        annotation.setIntergenic(true);
        if (!unspecStrandness) {    // strandness given, calculate distances on sense strand
            annotation.setgDist(calcDistancesNeighbours(firstOfPair, secondOfPair, it, start, end));
        } else {    // strandness not given, calculate distances on both strands
            annotation.setgDist(Math.min(calcDistancesNeighbours(firstOfPair, secondOfPair, it, start, end),
                    calcDistancesNeighbours(firstOfPair, secondOfPair, antiIt, start, end)));
        }
    }
}
