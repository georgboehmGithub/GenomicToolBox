import BamUtils.Annotation;
import GenomicAnnotations.Gene;
import GenomicAnnotations.GenomeAnnotations;
import GenomicAnnotations.Region;
import augmentedTree.IntervalTree;
import htsjdk.samtools.SAMRecord;

import java.util.ArrayList;
import java.util.HashMap;

public class BamRunner {
    private GenomeAnnotations gAnnotations;
    private HashMap<String, Integer> rawReadCounts;

    public BamRunner() {
        this.rawReadCounts = new HashMap<>();
    }

    public HashMap<String, Integer> getRawReadCounts() {
        return rawReadCounts;
    }

    public void setgAnnotations(GenomeAnnotations gAnnotations) {
        this.gAnnotations = gAnnotations;
    }
/*
    public static void main(String[] args) throws IOException {
        BamRunner runner = new BamRunner();
        Args argsParser = new Args();
        JCommander jct = JCommander.newBuilder()
                .addObject(argsParser)
                .build();
        jct.parse(args);
        if (argsParser.isHelp()) {
            jct.usage();
        } else {
            // frstrand
            String frStrand;
            boolean noFrStrand = false;
            if (!(argsParser.frstrand.equals("false") || argsParser.frstrand.equals("true"))) {
                noFrStrand = true;
            }

            // initializers
            GtfParser testParser = new GtfParser();
            GenomeAnnotations gAnnotations = testParser.parseGtf(new File(argsParser.inputGtf));
            runner.setgAnnotations(gAnnotations);
            File bamFile = new File(argsParser.bam);
            SamReader samReader = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
            Iterator<SAMRecord> bamIterator = samReader.iterator();

            // actual file work
            GeneralUtils BamGu = new GeneralUtils(gAnnotations);
            BamWriter bamWriter = new BamWriter(argsParser.outputPath);
            SAMRecord sr;
            String chromsome = "";
            boolean first = true;
            var lookup = new HashMap<String, SAMRecord>();
            String strandness;
            IntervalTree<Region> intervalTree = new IntervalTree<>();
            IntervalTree<Region> intervalTreeAnti = new IntervalTree<>();
            SAMRecord firstOfPair;
            SAMRecord secondOfPair;
            frStrand = noFrStrand ? "undefined" : argsParser.frstrand;
            int allReads = 0;

            while (bamIterator.hasNext()) {
                sr = bamIterator.next();
                // first SAMRecord
                if (!first) {
                    if (!chromsome.equals(sr.getReferenceName())) {
                        lookup.clear();
                        BamGu.getPcrCounts().clear();
                    }
                } else {
                    first = false;
                }
                String readName = sr.getReadName();
                if (BamGu.validRead(sr)) {
                    if (!lookup.containsKey(readName)) {
                        lookup.put(readName, sr);
                    } else {
                        SAMRecord mate = lookup.get(readName);

                        lookup.remove(readName);    // read pair irrelevant now

                        firstOfPair = sr.getFirstOfPairFlag() ? sr : mate;
                        secondOfPair = !sr.getFirstOfPairFlag() ? sr : mate;
                        strandness = runner.getStrandness(frStrand, firstOfPair, secondOfPair);
                        BamGu.setAnnotation(firstOfPair, secondOfPair, strandness);
                        intervalTree = runner.getIntervalTree(noFrStrand, sr, intervalTree, strandness);
                        intervalTreeAnti = runner.getIntervalTreeAnti(noFrStrand, sr, intervalTreeAnti, strandness);

                        if (BamGu.validReadPair(intervalTree, intervalTreeAnti, noFrStrand)) {
                            // further handle read pair except split inconsistency
                            if (!BamGu.splitConsistent(firstOfPair, secondOfPair)) {
                                // handle split inconsistency
                                if (argsParser.plotting.equals(""))
                                    bamWriter.writeSI(BamGu.getAnnotation());
                            } else {
                                // handle not split inconsistent read pair
                                BamGu.readpairMismatchCount(firstOfPair, secondOfPair);
                                BamGu.readPairClipping(firstOfPair, secondOfPair);
                                int pcrIndex = BamGu.calculatePcrIndices(firstOfPair, secondOfPair, strandness);

                                if (BamGu.getAnnotation().isIntergenic()) {
                                    if (argsParser.plotting.equals(""))
                                        bamWriter.writeIntergenic(BamGu.getAnnotation(), pcrIndex);
                                } else {
                                    allReads++;
                                    // plot data if necessary and flag given
                                    if (!argsParser.plotting.equals("")) {
                                        runner.plotData(BamGu.getAnnotation(), pcrIndex, false);
                                    } else {
                                        bamWriter.writeGenic(BamGu.getAnnotation(), pcrIndex);
                                    }
                                }
                            }
                        }
                    }
                }
                chromsome = sr.getReferenceName();  // reference = chromosome
                BamGu.setAnnotation(BamGu.getAnnotation().clear());  // holds all relevant information
            }
            if (!argsParser.plotting.equals(""))
                bamWriter.writeRPKMs(runner.getRawReadCounts(), gAnnotations.getGenes(), allReads);
            bamWriter.close();
        }
    }
*/
    // TODO

    /**
     * @param annotation
     * @param pcrIndex
     * @param allReads
     */
    public void plotData(Annotation annotation, int pcrIndex, boolean allReads) {
        if (!allReads && pcrIndex != 0) {
            return;
        }
        HashMap<Integer, ArrayList<Gene>> allGenes = annotation.getRelGenes();
        ArrayList<Gene> relGenes = allGenes.get(annotation.getHighestFeature());
        for (Gene curGene : relGenes) {
            String curId = curGene.getId();
            if (getRawReadCounts().containsKey(curId)) {
                int curCount = getRawReadCounts().get(curId);
                getRawReadCounts().put(curId, curCount + 1);
            } else {
                getRawReadCounts().put(curId, 1);
            }
        }
    }

    public String getStrandness(String frstrand, SAMRecord firstOfPair, SAMRecord secondOfPair) {
        if (frstrand.equals("true")) {
            return firstOfPair.getReadNegativeStrandFlag() ? "-" : "+";
        } else if (frstrand.equals("false")) {
            return secondOfPair.getReadNegativeStrandFlag() ? "-" : "+";
        } else {
            return "unspecific";
        }
    }

    public IntervalTree<Region> getIntervalTree(boolean noFrStrand, SAMRecord sr,
                                                IntervalTree<Region> intervalTree, String strandness) {
        if (!noFrStrand) {
            intervalTree = gAnnotations.getIntervalTrees().get(sr.getReferenceName() + strandness);
        } else {
            intervalTree = gAnnotations.getIntervalTrees().get(sr.getReferenceName() + "+");
        }
        return intervalTree;
    }

    public IntervalTree<Region> getIntervalTreeAnti(boolean noFrStrand, SAMRecord sr,
                                                    IntervalTree<Region> intervalTreeAnti, String strandness) {
        if (!noFrStrand) {
            intervalTreeAnti = strandness.equals("+") ?
                    gAnnotations.getIntervalTrees().get(sr.getReferenceName() + "-") :
                    gAnnotations.getIntervalTrees().get(sr.getReferenceName() + "+");
        } else {
            intervalTreeAnti = gAnnotations.getIntervalTrees().get(sr.getReferenceName() + "-");
        }
        return intervalTreeAnti;
    }
}
