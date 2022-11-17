package BamUtils;

import GenomicAnnotations.Gene;
import GenomicAnnotations.Region;
import GenomicAnnotations.RegionVector;
import GenomicAnnotations.Transcript;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class BamWriter {
    private BufferedWriter bw;

    public BamWriter(String outputPath) throws IOException {
        this.bw = new BufferedWriter(new FileWriter(outputPath));
    }

    public void writeRPKMs(HashMap<String, Integer> rawReadCounts, HashMap<String, Gene> allGenes, int allReads) throws IOException {
        var RPKMS = new HashMap<Double, Integer>();
        ArrayList<Region> allTranscriptsExons;

        for (Map.Entry<String, Integer> entry : rawReadCounts.entrySet()) { // counts for individual genes
            allTranscriptsExons = new ArrayList<>();
            Gene curGene = allGenes.get(entry.getKey());
            int readCount = entry.getValue();
            // current genes transcript set
            for (Map.Entry<String, Transcript> transcriptEntry : curGene.getTranscripts().entrySet()) {
                Transcript curTranscript = transcriptEntry.getValue();
                allTranscriptsExons.addAll(curTranscript.getExons().getRegions());
            }
            RegionVector mergedTranscripts = new RegionVector(allTranscriptsExons);
            mergedTranscripts = mergedTranscripts.mergedRv();

            double transcriptsLength = 0.0;
            for (Region reg : mergedTranscripts.getRegions()) {
                transcriptsLength += reg.getLength();
            }
            double RPKM = (readCount / (allReads / 1000000.0)) / (transcriptsLength / 1000);
                if (RPKMS.containsKey(RPKM)) {
                    int curRPKMCount = RPKMS.get(RPKM);
                    RPKMS.put(RPKM, curRPKMCount + 1);
                } else {
                    RPKMS.put(RPKM, 1);
                }
        }
        bw.write("RPKM" + "\t" + "reads" + "\n");
        for (Map.Entry<Double, Integer> RPKM : RPKMS.entrySet()) {
            bw.write(RPKM.getKey() + "\t" + "non_pcr_duplicated_reads" + "\n");
        }
    }

    public void writeSI(Annotation annot) throws IOException {
        bw.write(annot.toStringSplitInconsistent());
    }

    public void writeIntergenic(Annotation annot, int pcrIndex) throws IOException {
        bw.write(annot.toStringIntergenic(pcrIndex));
    }

    // TODO in Annotation Klasse toString methoden verlagern alles
    public void writeGenic(Annotation annot, int pcrIndex) throws IOException {
        HashMap<Integer, ArrayList<Gene>> allGenes = annot.getRelGenes();
        ArrayList<Gene> relGenes;
        HashMap<String, ArrayList<Transcript>> relTranscripts;
        relGenes = allGenes.get(annot.getHighestFeature());
        TreeSet<String> geneIds = new TreeSet<>();

        for (Gene gene : relGenes) {
            geneIds.add(gene.getId());
        }
        int gcount = geneIds.size();
        TreeSet<String> uniqTransIds = new TreeSet<>();

        String matchedGenes = "";
        if (annot.getHighestFeature() == 2) {
            relTranscripts = annot.getRelTranscripts();
            for (Map.Entry<String, ArrayList<Transcript>> entry : relTranscripts.entrySet()) {
                for (Transcript t : entry.getValue()) {
                    uniqTransIds.add(t.getId());
                }
            }
            gcount = 0; // potentially not the same as geneIds.size()
            for (Gene curGene : relGenes) {
                if (geneIds.contains(curGene.getId())) {
                    geneIds.remove(curGene.getId());
                    String geneEntry = "";
                    geneEntry += curGene.getId() + "," + curGene.getBiotype() + ":";
                    if (uniqTransIds.size() > 0) {
                        gcount++;
                        for (Transcript curTrans : relTranscripts.get(curGene.getId())) {
                            if (uniqTransIds.contains(curTrans.getId())) {
                                uniqTransIds.remove(curTrans.getId());
                                geneEntry += curTrans.getId() + ",";
                            }
                        }
                        geneEntry = geneEntry.substring(0, geneEntry.length() - 1);
                        matchedGenes += geneEntry;
                        if (geneIds.size() > 0) {
                            matchedGenes += "|";
                        }
                    }
                }
            }
            if (matchedGenes.charAt(matchedGenes.length() - 1) == '|') {
                matchedGenes = matchedGenes.substring(0, matchedGenes.length() - 1);
            }
        } else if (annot.getHighestFeature() == 1) {
            for (Gene curGene : relGenes) {
                if (geneIds.contains(curGene.getId())) {
                    geneIds.remove(curGene.getId());
                    matchedGenes += curGene.getId() + "," + curGene.getBiotype() + ":" + "MERGED";
                    if (geneIds.size() > 0) {
                        matchedGenes += "|";
                    }
                }
            }
        } else {    // highest feature == 0
            for (Gene curGene : relGenes) {
                if (geneIds.contains(curGene.getId())) {
                    geneIds.remove(curGene.getId());
                    matchedGenes += curGene.getId() + "," + curGene.getBiotype() + ":" + "INTRON";
                    if (geneIds.size() > 0) {
                        matchedGenes += "|";
                    }
                }
            }
        }
        bw.write(annot.getFirstOfPair().getReadName() + "\t" + "mm:" +
                annot.getMismatches() + "\t" + "clipping:" + annot.getClipped() +
                "\t" + "nsplit:" + annot.getSplitCount() + "\t" + "gcount:" + gcount + "\t"
                + matchedGenes + "\t" + "pcrindex:" + pcrIndex + "\n");
    }

    public void close() throws IOException {
        this.bw.close();
    }
}
