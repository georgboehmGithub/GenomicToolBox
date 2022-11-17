import Events.ExonSkippingEvent;
import GenomicAnnotations.Gene;
import GenomicAnnotations.GenomeAnnotations;
import augmentedTree.IntervalTree;

import java.io.*;
import java.util.HashMap;

public class GtfParser {

    /**
     * Iterates over input file line by line and parses each entry into its
     * corresponding "feature" object.
     *
     * @param gtf        input file in gtf format
     * @return GenomeAnnotations object containing all parsed objects
     * @throws IOException
     */
    public GenomeAnnotations parseGtf(File gtf) throws IOException{
        GenomeAnnotations GtfAnnotations = new GenomeAnnotations();
        try (BufferedReader br = new BufferedReader(new FileReader(gtf))) {
            String line;

            while ((line = br.readLine()) != null) {
                if (!line.startsWith("#")) {
                    String[] entry = splitString(line, '\t');
                    HashMap<String, String> attributes = generateAttributesHM(entry);

                    String preGeneId = attributes.get("gene_id");
                    String geneId = preGeneId.substring(1, preGeneId.length() - 1);
                    String feature = entry[2];

                    // interval trees
                    if (!GtfAnnotations.getIntervalTrees().containsKey(entry[0] + GtfAnnotations.identifyStrand(entry))) {
                        GtfAnnotations.getIntervalTrees().put(entry[0] + GtfAnnotations.identifyStrand(entry), new IntervalTree<>());
                    }

                    switch (feature) {
                        case "gene":
                            if (!GtfAnnotations.getGenes().containsKey(geneId)) {
                                GtfAnnotations.parseGene(entry, geneId, attributes);
                            }
                            break;
                        case "transcript":
                            GtfAnnotations.parseTranscript(entry, geneId, attributes);
                            break;
                        case "exon": {
                            String preTranscriptId = attributes.get("transcript_id");
                            String transcriptId = preTranscriptId.substring(1, preTranscriptId.length() - 1);
                            GtfAnnotations.parseExon(entry, geneId, transcriptId, attributes);
                            break;
                        }
                        case "CDS": {
                            String preTranscriptId = attributes.get("transcript_id");
                            String transcriptId = preTranscriptId.substring(1, preTranscriptId.length() - 1);
                            GtfAnnotations.parseCds(entry, geneId, transcriptId, attributes);
                            break;
                        }
                    }
                }
            }
        } catch (IOException e) {
            throw new IOException("Error while reading gtf input file");
        }
        //GtfAnnotations.getExonSkippings();
        //if (outputPath.length() > 0) writeEsFile(GtfAnnotations, outputPath);
        return GtfAnnotations;
    }

    public String[] splitString(String line, final char delimiter) {
        CharSequence[] temp = new CharSequence[(line.length() / 2) + 1];
        int wordCount = 0;
        int i = 0;
        int j = line.indexOf(delimiter, 0); // first substring

        while (j >= 0) {
            temp[wordCount++] = line.substring(i, j);
            i = j + 1;
            j = line.indexOf(delimiter, i); // rest of substrings
        }

        temp[wordCount++] = line.substring(i); // last substring

        String[] result = new String[wordCount];
        System.arraycopy(temp, 0, result, 0, wordCount);

        return result;
    }

    /**
     * Reads the last attributes column of a gtf file into a Hashmap
     * containing attribute keys with corresponding values
     *
     * @param entry line being read
     * @return attributes hashmap for further use
     */
    public HashMap<String, String> generateAttributesHM(String[] entry) {
        HashMap<String, String> attributes = new HashMap<>();
        String[] aStrings = splitString(entry[8], ';');

        for (String attribute : aStrings) {
            String[] attKeyVal = splitString(attribute, ' ');
            if (attKeyVal.length == 2) attributes.put(attKeyVal[0], attKeyVal[1]);
            else if (attKeyVal.length == 3) attributes.put(attKeyVal[1], attKeyVal[2]);
        }
        return attributes;
    }

    /**
     * Writes exon skipping events into output file, line by line.
     *
     * @param annotations entirety of the gtf parsed into objects
     * @param outputPath  outputfile path
     */
    public void writeEsFile(GenomeAnnotations annotations, String outputPath) {
        try {
            FileWriter fos = new FileWriter(outputPath);
            PrintWriter dos = new PrintWriter(fos);
            dos.println("id\tsymbol\tchr\tstrand\tnprots\tntrans\tSV\tWT\tWT_prots\tSV_prots\tmin_skipped_exon" +
                    "\tmax_skipped_exon\tmin_skipped_bases\tmax_skipped_bases");
            for (ExonSkippingEvent esEvent : annotations.getExonSkippings()) {
                Gene esGene = esEvent.getGene();
                dos.print(esGene.toString());
                dos.print(esEvent.toString());
                dos.println();
            }
            dos.close();
            fos.close();
        } catch (IOException e) {
            System.out.println("Error Printing Tab Delimited File");
        }
    }
/*
    public static void main(String[] args) throws IOException {
        Args parser = new Args();
        JCommander jct = JCommander.newBuilder()
                .addObject(parser)
                .build();
        jct.parse(args);
        if (parser.isHelp()) {
            jct.usage();
        } else {
            String filePath = parser.inputGtf;
            String outputPath = parser.outputPath;
            if (filePath.equals("") || outputPath.equals("")) {
                jct.usage();
            } else {
                File sampleFile = new File(filePath);
                GtfParser testParser = new GtfParser();
                //GenomeAnnotations annotations = testParser.parseGtf(sampleFile, outputPath);
            }
        }
    }
 */
}
