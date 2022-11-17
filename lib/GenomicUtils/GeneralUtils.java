package GenomicUtils;

import GenomicAnnotations.Gene;
import GenomicAnnotations.Region;
import GenomicAnnotations.RegionVector;
import GenomicAnnotations.Transcript;

import java.io.*;
import java.util.*;
import java.util.concurrent.ThreadLocalRandom;

public class GeneralUtils {
    public HashMap<Character, char[]> mutationAccess;
    public double mutationRate;
    private Random rand;

    public GeneralUtils() {
        fillMutationAccess();
        this.rand = new Random();
    }

    public GeneralUtils(String mutationRate) {
        fillMutationAccess();
        this.mutationRate = calcMutationRate(mutationRate);
        this.rand = new Random();
    }

    public void fillMutationAccess() {
        this.mutationAccess = new HashMap<>();
        this.mutationAccess.put('C', new char[]{'A', 'T', 'G'});
        this.mutationAccess.put('A', new char[]{'C', 'T', 'G'});
        this.mutationAccess.put('G', new char[]{'A', 'T', 'C'});
        this.mutationAccess.put('T', new char[]{'A', 'C', 'G'});
    }

    /**
     * Example: "1.0%" input, 0.01 output.
     * @param mutationRate
     * @return
     */
    public double calcMutationRate(String mutationRate) {
        mutationRate = mutationRate.substring(0, mutationRate.length() - 1);
        double parsedMr = ((int) Double.parseDouble(mutationRate)) / 100.0;
        return parsedMr;
    }

    /**
     * Retrieves transcript sequence of given gene sequence by
     * cutting the exons out of the sequencec and concatenating them.
     * @param gene sequence origin
     * @param transcript transcript under observation
     * @param gse does the extraction work
     * @return transcript sequence, already respecting strandness
     */
    public String getTransSeq(Gene gene, Transcript transcript, GenomeSequenceExtractor gse) {
        String geneSeq = gse.getSequence(gene.getChromosome(), gene.getStart(), gene.getEnd());
        String transSeq = "";
        transcript.getExons().sortRegionVector();
        for (Region exon : transcript.getExons().getRegions()) {
            int seqStart = exon.getStart() - gene.getStart();
            int seqEnd = seqStart + (exon.getEnd() - exon.getStart() + 1);
            transSeq += geneSeq.substring(seqStart, seqEnd);
        }
        if (gene.getStrand().equals("-")) {
            transSeq = revComp(transSeq);
        }
        return transSeq;
    }

    public RegionVector getCoveredRegions(Transcript trans, FragmentData fd, boolean fw, int transLength, String curStrand) {
        RegionVector gCoveredRegions = new RegionVector();
        ArrayList<Region> transExons = trans.getExons().getRegions();
        Region curExon = transExons.get(0);
        int gStart;
        int gStop;

        gStart = fw ? fd.getfLocalRegion().getStart() : fd.getfLocalRegion().getEnd() - fd.getReadLength();
        gStop = gStart + fd.getReadLength();

        if (curStrand.equals("-")) {
            gStart = transLength - gStop;
        }

        int stretch = 0;
        int i = 0;

        //find starting exon
        while (i < transExons.size()) {
            curExon = transExons.get(i);
            if (stretch + curExon.getLength() < gStart + 1) {
                i++;
                stretch += curExon.getLength();
            } else {
                break;
            }
        }
        //find other exons
        int startPos = curExon.getStart() + (gStart - stretch);
        int basesLeft = fd.getReadLength();
        //i = i == transExons.size() ? i - 1 : i;

        while (basesLeft > 0) {
            if (basesLeft < curExon.getEnd() - startPos + 1) {
                gCoveredRegions.getRegions().add(new Region(startPos, startPos + basesLeft));
                break;
            } else {
                gCoveredRegions.getRegions().add(new Region(startPos, curExon.getEnd() + 1));
                basesLeft -= (curExon.getEnd() + 1 - startPos);
            }
            i++;
            if (i == transExons.size()) break;
            curExon = transExons.get(i);
            startPos = curExon.getStart();
        }
        return gCoveredRegions;
    }

    public RegionVector getCoveredRegions2(Transcript trans, FragmentData fd, boolean fw) {
        RegionVector gCoveredRegions = new RegionVector();
        ArrayList<Region> transExons = trans.getExons().getRegions();
        Region curExon;
        Region coveredRegion;

        int readStart = fw ? fd.getfLocalRegion().getStart() + 1 : fd.getfLocalRegion().getEnd() - fd.getReadLength();
        int basesLeft = fd.getReadLength();
        int stretch = transExons.get(0).getLength();
        int lastStretch = 0;

        for (int i = 0; i < transExons.size(); i++) {
            if (basesLeft <= 0) {
                break;
            }
            curExon = transExons.get(i);
            if (i > 0) {
                stretch += (curExon.getLength());
            }
            if (basesLeft < fd.getReadLength()) {
                if (curExon.getStart() + basesLeft <= curExon.getEnd()) {
                    coveredRegion = new Region(curExon.getStart(), curExon.getStart() + basesLeft);
                    gCoveredRegions.getRegions().add(coveredRegion);
                    break;
                } else {
                    coveredRegion = new Region(curExon.getStart(), curExon.getEnd());
                    gCoveredRegions.getRegions().add(coveredRegion);
                    basesLeft -= coveredRegion.getLength() - 1;
                }
            } else {
                if (stretch > readStart) {
                    int startPos = readStart - lastStretch; //TODO: -1??
                    if (curExon.getStart() + startPos + basesLeft <= curExon.getEnd()) {
                        coveredRegion = new Region(curExon.getStart() + startPos, curExon.getStart() + startPos + basesLeft);
                        gCoveredRegions.getRegions().add(coveredRegion);
                        break;
                    } else {
                        coveredRegion = new Region(curExon.getStart() + startPos, curExon.getEnd());
                        gCoveredRegions.getRegions().add(coveredRegion);
                        basesLeft -= coveredRegion.getLength() - 1;
                    }
                } else if (stretch == readStart) { // Only 1 overlap
                    coveredRegion = new Region(curExon.getEnd(), curExon.getEnd()); //TODO look over this
                    gCoveredRegions.getRegions().add(coveredRegion);
                    basesLeft -= 1;
                }
            }
            lastStretch = stretch;
        }
        return gCoveredRegions;
    }

    public FragmentData getFragmentData(int mean, int sd, int readLength, String transcriptSeq,
                                        FragmentData fragmentData, String curStrand) {
        String[] fwMutData;
        String[] rvMutData;
        int transLength = transcriptSeq.length();
        if (transLength < readLength) {
            fragmentData.setInvalidSequence(true);
            return fragmentData;
        }

        int ndVal = (int) Math.round(rand.nextGaussian() * sd + mean);
        while (ndVal < readLength) {
            ndVal = (int) Math.round(rand.nextGaussian() * sd + mean);
        }
        int fragLength = ndVal;

        if (fragLength >= transLength) {
            fragLength = transLength - 1;
        }

        int fStart = ThreadLocalRandom.current().nextInt(transLength - fragLength);

        String fragment = transcriptSeq.substring(fStart, fStart + fragLength); //TODO: end is exclusive!!! Keep in mind

        fwMutData = mutated(fragment.substring(0, readLength));
        rvMutData = mutated(revComp(fragment.substring(fragment.length() - readLength, fragment.length())));

        fragmentData.setReads(fwMutData[0], fwMutData[1], rvMutData[0], rvMutData[1]);
        fragmentData.setfLocalRegion(new Region(fStart, fStart + fragLength));
        return fragmentData;
    }

    /**
     * Builds the reverse complement of given sequence.
     * @param sequence sequence to rev. comp.
     * @return rev. comp.
     */
    public String revComp(String sequence) {
        var rc = new StringBuilder();
        char c;
        for (int i = sequence.length() - 1; i >= 0; i--) {
            c = sequence.charAt(i);
            switch (c) {
                case 'A':
                    rc.append('T');
                    break;
                case 'T':
                    rc.append('A');
                    break;
                case 'C':
                    rc.append('G');
                    break;
                case 'G':
                    rc.append('C');
                    break;
            }
        }
        return rc.toString();
    }

    /**
     * Simulates pseudo random mutations of bases for a genomic sequence.
     *
     * @param template string to mutate
     * @return mutated template string
     */
    public String[] mutated(String template) {
        String mutations = "";
        char c;
        int probability;
        int randomAccess;
        String mutatedString = "";
        int probabilityFactor = (int) (1.0 / mutationRate);
        int luckyNumber = (int) (Math.random() * probabilityFactor);

        for (int i = 0; i < template.getBytes().length; i++) {
            c = template.charAt(i);
            probability = (int) (Math.random() * probabilityFactor);
            if (probability == luckyNumber) {
                randomAccess = (int) (Math.random() * 2);
                mutatedString += mutationAccess.get(c)[randomAccess];
                mutations += i + ",";
            } else {
                mutatedString += c;
            }
        }
        if (mutations.length() > 0) {
            mutations = mutations.substring(0, mutations.length() - 1);
        }
        return new String[]{mutatedString, mutations};
    }

    public void barplotData(String rmInfo) {
        try (BufferedReader br = new BufferedReader(new FileReader(rmInfo))) {
            String line;
            String[] data;
            String[] fwRegions;
            String[] rwRegions;
            boolean mismatches = false;
            int nonSplitCount = 0;
            int nonSplitWoMm = 0;
            int splitCount = 0;
            int splitWoMm = 0;
            boolean atleastFive = true;
            int splitWoMmFive = 0;
            int readCount = 0;

            boolean header = true;
            while ((line = br.readLine()) != null) {
                if (header) {
                    header = false;
                } else {
                    readCount++;
                    data = line.split("\t");
                    if (data.length >= 9) mismatches = true;
                    fwRegions = data[6].split("\\|");
                    rwRegions = data[7].split("\\|");

                    if (fwRegions.length > 1 || rwRegions.length > 1) { // split read-pair
                        splitCount++;
                        if (!mismatches) {
                            splitWoMm++;
                            String[] regionSa;
                            for (String regionS : fwRegions) {
                                regionSa = regionS.split("-");
                                // Check wether length of each region is atleast 5 long
                                Long end = Long.parseLong(regionSa[1]);
                                Long start = Long.parseLong(regionSa[0]);
                                if (end - start + 1 < 5) {
                                    atleastFive = false;
                                    break;
                                }
                            }
                            if (atleastFive) {  // If it still hold true for fw, check for rw
                                for (String regionS : rwRegions) {
                                    regionSa = regionS.split("-");
                                    // Check wether length of each region is atleast 5 long
                                    Long end = Long.parseLong(regionSa[1]);
                                    Long start = Long.parseLong(regionSa[0]);
                                    if (end - start + 1 < 5) {
                                        atleastFive = false;
                                        break;
                                    }
                                }
                            }
                            if (atleastFive) splitWoMmFive++;
                        }
                    } else {    // non-split read-pair
                        nonSplitCount++;
                        if (!mismatches) nonSplitWoMm++;
                    }
                    mismatches = false;
                    atleastFive = true;
                }
            }
            int[] plotData = new int[]{readCount, nonSplitCount, nonSplitWoMm, splitCount, splitWoMm,
                    splitWoMmFive};
            String[] ids = new String[]{"allReads", "nonSplitCount", "nonSplitWoMm", "splitCount", "splitWoMm",
                    "splitWoMmFive"};
            System.out.println();
            BufferedWriter bwBarplots = new BufferedWriter(new FileWriter("readsData.csv"));
            bwBarplots.write("objective" + "\t" + "occurences" + "\n");
            for (int i = 0; i <= 5; i++) {
                bwBarplots.write(ids[i] + "\t" + plotData[i] + "\n");
            }
            bwBarplots.close();
            } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void plotDataFL(String rmInfo) {
        try (BufferedReader br = new BufferedReader(new FileReader(rmInfo))) {
            String line;
            String[] data;
            String tFwRegvec;
            String tRwRegvec;
            boolean header = true;
            int fStart;
            int fEnd;
            int fLength;
            var flDistribution = new TreeMap<Integer, Integer>();
            while ((line = br.readLine()) != null) {
                if (header) {
                    header = false;
                } else {
                    data = line.split("\t");

                    tFwRegvec = data[4];
                    tRwRegvec = data[5];

                    fStart = Integer.parseInt(tFwRegvec.split("-")[0]);
                    fEnd = Integer.parseInt(tRwRegvec.split("-")[1]);
                    fLength = fEnd - fStart;
                    if (flDistribution.containsKey(fLength)) {
                        int val = flDistribution.get(fLength);
                        flDistribution.put(fLength, val + 1);
                    } else {
                        flDistribution.put(fLength, 1);
                    }
                }
            }
            BufferedWriter bwFl = new BufferedWriter(new FileWriter("flDistribution.csv"));
            bwFl.write("flength" + "\t" + "counts" + "\n");
            for (Map.Entry<Integer,Integer> entry : flDistribution.entrySet()) {
                int flength = entry.getKey();
                int counts = entry.getValue();
                bwFl.write(flength + "\t" + counts + "\n");
            }
            bwFl.close();
            } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void plotDataMt(String rmInfo) {
        try (BufferedReader br = new BufferedReader(new FileReader(rmInfo))) {
            String line;
            String[] data;
            String[] fwMut;
            String[] rwMut;
            int mutationCount;
            int occurences;

            boolean header = true;
            var mutDistribution = new TreeMap<Integer, Integer>();
            while ((line = br.readLine()) != null) {
                if (header) {
                    header = false;
                } else {
                    data = line.split("\t");
                    if (data.length == 10) {
                        fwMut = data[8].split(",");
                        rwMut = data[9].split(",");
                        if (fwMut[0].equals("")) {
                            mutationCount = rwMut.length;
                        } else {
                            mutationCount = fwMut.length + rwMut.length;
                        }
                    } else if (data.length == 9) {
                        fwMut = data[8].split(",");
                        mutationCount = fwMut.length;
                    } else {
                        mutationCount = 0;
                    }
                    if (mutDistribution.containsKey(mutationCount)) {
                        occurences = mutDistribution.get(mutationCount);
                        mutDistribution.put(mutationCount, occurences + 1);
                    } else {
                        mutDistribution.put(mutationCount, 1);
                    }
                }
            }
            BufferedWriter bwMut = new BufferedWriter(new FileWriter("mutDistribution.csv"));
            bwMut.write("mutationCount" + "\t" + "occurences" + "\n");
            for (Map.Entry<Integer,Integer> entry : mutDistribution.entrySet()) {
                bwMut.write(entry.getKey() + "\t" + entry.getValue() + "\n");
            }
            bwMut.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        GeneralUtils g = new GeneralUtils();
        String mappingInfo = "C:\\Users\\georg\\IdeaProjects\\assignment1\\src\\read.mappinginfo";
        // fragment lengths
        // g.plotDataFL(mappingInfo);

        // mutation distribution
        //g.plotDataMt(mappingInfo);

        // barplot data
        g.barplotData(mappingInfo);
    }
}
