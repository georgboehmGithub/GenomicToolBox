package GenomicAnnotations;

import Events.ExonSkippingEvent;
import augmentedTree.IntervalTree;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class GenomeAnnotations {
    private final HashMap<String, Gene> genes;
    private final HashMap<String, Transcript> transcripts;
    private final ArrayList<ExonSkippingEvent> esEvents;
    private final HashMap<String, IntervalTree<Region>> intervalTrees;

    public GenomeAnnotations() {
        this.genes = new HashMap<>();
        this.transcripts = new HashMap<>();
        this.esEvents = new ArrayList<>();
        this.intervalTrees = new HashMap<>();
    }

    /**
     * Iterates over each column relevant for identifying a entries strand.
     *
     * @param entry line being read
     * @return corresponding strand for a given entry
     */
    public String identifyStrand(String[] entry) {
        if (entry[5].equals("-") | entry[6].equals("-") | entry[7].equals("-")) {
            return "-";
        } else if (entry[5].equals("+") | entry[6].equals("+") | entry[7].equals("+")) {
            return "+";
        } else {
            return "undefined";
        }
    }

    /**
     * Parses either a line with feature "gene" into a gene object, or
     * creates a not yet initialized gene for a different feature line (e.g. exon, transcript)
     *
     * @param entry      current line being read
     * @param geneId     the new gene's id
     * @param attributes the new gene's attributes
     */
    public void parseGene(String[] entry, String geneId, HashMap<String, String> attributes) {
        String chromosome = entry[0];
        String strand = identifyStrand(entry);
        String geneName = "";
        String biotype = "";


        if (attributes.containsKey("gene_name")) {
            String preGeneName = attributes.get("gene_name");
            geneName = preGeneName.substring(1, preGeneName.length() - 1);
        }
        if (attributes.containsKey("gene_biotype")) {
            String preBiotype = attributes.get("gene_biotype");
            biotype = preBiotype.substring(1, preBiotype.length() - 1);
        }
        if (geneName.length() == 0) geneName = "undefined";
        if (biotype.length() == 0) biotype = "undefined";

        int geneStart = Integer.parseInt(entry[3]);
        int geneEnd = Integer.parseInt(entry[4]);

        intervalTrees.get(chromosome + strand).add(new Region(geneStart, geneEnd, geneId));

        Gene newGene = new Gene(chromosome, geneId, geneName, biotype, strand,
                geneStart, geneEnd);
        this.getGenes().put(geneId, newGene);
    }

    /**
     * Parses either a line with feature "transcript" into a transcript object, or
     * creates a not yet initialized transcript for a different feature line (e.g. exon, cds)
     *
     * @param entry      current line being read
     * @param geneId     the current gene's id
     * @param attributes the new transcript's attributes
     */
    public void parseTranscript(String[] entry, String geneId, HashMap<String, String> attributes) {
        if (!genes.containsKey(geneId)) {
            parseGene(entry, geneId, attributes);
        }
        Gene currentGene = genes.get(geneId);
        String transcriptId = "";
        String transcriptName = "";

        if (attributes.containsKey("transcript_id")) {
            String preTranscriptId = attributes.get("transcript_id");
            transcriptId = preTranscriptId.substring(1, preTranscriptId.length() - 1);
        }
        if (attributes.containsKey("transcript_name")) {
            String preTranscriptName = attributes.get("transcript_name");
            transcriptName = preTranscriptName.substring(1, preTranscriptName.length() - 1);
        }
        if (transcriptId.length() == 0) transcriptId = "undefined";
        if (transcriptName.length() == 0) transcriptName = "undefined";

        int start = Integer.parseInt(entry[3]);
        int end = Integer.parseInt(entry[4]);
        Transcript newTranscript = new Transcript(transcriptId, transcriptName, start, end);
        currentGene.getTranscripts().put(transcriptId, newTranscript);
        transcripts.put(transcriptId, newTranscript);
    }

    /**
     * Parses a line with feature "exon" into a new Region, adding it to its corresponding transcript.
     *
     * @param entry        current line being read
     * @param geneId       the current gene's id
     * @param transcriptId the current transcript's id
     * @param attributes   the new transcript's attributes, if need to be created
     */
    public void parseExon(String[] entry, String geneId, String transcriptId, HashMap<String, String> attributes) {
        if (!transcripts.containsKey(transcriptId)) {
            parseTranscript(entry, geneId, attributes);
        }
        Region newExonRegion = new Region(Integer.parseInt(entry[3]), Integer.parseInt(entry[4]));
        Gene curGene = genes.get(geneId);
        curGene.getTranscripts().get(transcriptId).addExon(newExonRegion);
    }

    /**
     * Parses a line with feature "cds" into a new CDS Object, adding it to its corresponding transcript.
     *
     * @param entry        current line being read
     * @param geneId       the current gene's id
     * @param transcriptId the current transcript's id
     * @param attributes   the new cds attributes
     */
    public void parseCds(String[] entry, String geneId, String transcriptId, HashMap<String, String> attributes) {
        Region newCdsRegion = new Region(Integer.parseInt(entry[3]), Integer.parseInt(entry[4]));
        String protId = "undefined";
        if (!transcripts.containsKey(transcriptId)) {
            parseTranscript(entry, geneId, attributes);
        }
        if (attributes.containsKey("protein_id")) {
            String preProteinId = attributes.get("protein_id");
            protId = preProteinId.substring(1, preProteinId.length() - 1);
        }
        Gene curGene = genes.get(geneId);
        Transcript curTranscript = curGene.getTranscripts().get(transcriptId);

        if (curTranscript.getCDS() == null) {
            CDS newCDS = new CDS(protId, curTranscript.getId());
            newCDS.addCds(newCdsRegion);

            // Set vals for transcript
            curTranscript.setCodingSequence(newCDS);
            curTranscript.setProtId(protId);

            // Set vals for gene
            curGene.setNprots(curGene.getNprots() + 1);
            curGene.getCDSs().put(protId, newCDS);
        } else {
            curTranscript.getCDS().addCds(newCdsRegion);
        }
    }

    /**
     * Looks for SV introns in each gene of the gtf file. If an exon skipping event took place, a new
     * exonSkippingEvent object is initialized and calculated.
     *
     * @return GenomeAnnotation object's exon Skipping Events Arraylist if already calculates. Otherwise
     * calculates it and returns it.
     */
    public ArrayList<ExonSkippingEvent> getExonSkippings() {
        if (esEvents.size() == 0) {
            for (Map.Entry<String, Gene> geneEntry : genes.entrySet()) {
                Gene currentGene = geneEntry.getValue();
                HashSet<Region> candIntrons = currentGene.uniqueIntrons();
                for (Region candSVIntron : candIntrons) {
                    HashSet<String> SV = new HashSet<>();
                    HashSet<String> wtStart = new HashSet<>();
                    HashSet<String> wtEnd = new HashSet<>();
                    for (Map.Entry<String, Transcript> transcriptEntry : currentGene.getTranscripts().entrySet()) {
                        Transcript t = transcriptEntry.getValue();
                        String tId = t.getId();
                        ArrayList<Region> introns = t.getIntrons().getRegions();
                        for (Region intron : introns) {
                            if (intron.equals(candSVIntron)) {
                                SV.add(tId);
                            } else {
                                if (intron.getStart() == candSVIntron.getStart()) {
                                    wtStart.add(tId);
                                } else if (intron.getEnd() == candSVIntron.getEnd()) {
                                    wtEnd.add(tId);
                                }
                            }
                        }
                    }
                    HashSet<String> WT = new HashSet<>(wtStart);
                    WT.retainAll(wtEnd);
                    WT.removeAll(SV);
                    if (WT.size() > 0) {
                        ExonSkippingEvent newEsEvent = new ExonSkippingEvent(currentGene, SV, WT, candSVIntron);
                        newEsEvent.minMaxSkippedExons(WT, candSVIntron);
                        esEvents.add(newEsEvent);
                    }
                }
            }
        }
        return esEvents;
    }

    public void writePlotData(String outputPath, String fileName) {
        var uniqMaxExonCounts = new HashSet<Integer>();
        var data = new HashMap<Integer, ArrayList<ExonSkippingEvent>>();
        for (ExonSkippingEvent es : getExonSkippings()) {
            int posUniqCount = es.getMaxSkippedExons().size();
            uniqMaxExonCounts.add(posUniqCount);
        }
        for (int uniqCount : uniqMaxExonCounts) {
            data.put(uniqCount, new ArrayList<>());
        }
        for (ExonSkippingEvent es : getExonSkippings()) {
            int maxExonCount = es.getMaxSkippedExons().size();
            data.get(maxExonCount).add(es);
        }
        try {
            FileWriter fos = new FileWriter(outputPath);
            PrintWriter dos = new PrintWriter(fos);
            dos.println("maxSkippedExons\tesEvCount\tfilename");
            for (Map.Entry<Integer, ArrayList<ExonSkippingEvent>> dataPoint : data.entrySet()) {
                dos.print(dataPoint.getKey() + "\t");
                dos.print(dataPoint.getValue().size() + "\t");
                dos.print(fileName);
                dos.println();
            }
            dos.close();
            fos.close();
        } catch (IOException e) {
            System.out.println("Error Printing Tab Delimited File");
        }
    }

    public HashMap<String, IntervalTree<Region>> getIntervalTrees() {
        return intervalTrees;
    }

    public Transcript getTranscriptById(String id) {
        return transcripts.get(id);
    }

    public Gene getGeneById(String id) {
        return genes.get(id);
    }

    public HashMap<String, Gene> getGenes() {
        return genes;
    }
}
