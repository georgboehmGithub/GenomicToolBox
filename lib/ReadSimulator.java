import GenomicUtils.GeneralUtils;
import GenomicUtils.IndexFile;

import java.io.File;

public class ReadSimulator {
    private final File gtf;
    private final IndexFile fidx;
    private final File fastaGenome;
    private final File readcounts;
    private final String od;
    private final GeneralUtils gu;
    private final int frlength;
    private final int sd;
    private final int readlength;

    public ReadSimulator(String gtf, String fidx, String fg, String rc, String od,
                         String mutationrate, String frlength, String sd, String readlength) {
        this.gtf = new File(gtf);
        this.fidx = new IndexFile(fidx);
        this.fastaGenome = new File(fg);
        this.readcounts = new File(rc);
        this.od = od;
        this.gu = new GeneralUtils(mutationrate);
        this.frlength = Integer.parseInt(frlength);
        this.sd = Integer.parseInt(sd);
        this.readlength = Integer.parseInt(readlength);
    }
/*
    public static void main(String[] args) throws IOException {
        Args argsParser = new Args();
        JCommander jct = JCommander.newBuilder()
                .addObject(argsParser)
                .build();
        jct.parse(args);
        if (argsParser.isHelp()) {
            jct.usage();
        }

        ReadSimulator rs = new ReadSimulator(argsParser.inputGtf, argsParser.fidx,
                argsParser.fasta, argsParser.readcounts, argsParser.od, argsParser.mutationrate,
                argsParser.frlength, argsParser.SD, argsParser.length);

        GtfParser gtfP = new GtfParser();
        GenomeAnnotations annotations = gtfP.parseGtf(rs.gtf);
        GenomeSequenceExtractor gse = new GenomeSequenceExtractor(rs.fastaGenome, rs.fidx);

        try (BufferedReader br = new BufferedReader(new FileReader(rs.readcounts))) {
            boolean header = true;
            int readcount;
            String line;
            String[] entry;
            String curTransSeq;
            String curStrand;
            String qualityScore = "I".repeat(rs.readlength);
            FragmentData fragmentData = new FragmentData(rs.readlength);

            BufferedWriter bwMappingInfo = new BufferedWriter(new FileWriter(rs.od + "/" + "read.mappinginfo"));
            bwMappingInfo.write("readid" + "\t" + "chr" + "\t" + "gene" + "\t" + "transcript" + "\t" + "t_fw_regvec" + "\t"
                    + "t_rw_regvec" + "\t" + "fw_regvec" + "\t" + "rw_regvec" + "\t" + "fw_mut" + "\t" + "rw_mut" + "\n");

            MappingWriter mw = new MappingWriter(fragmentData, bwMappingInfo);
            BufferedWriter fwFastqWriter = new BufferedWriter(new FileWriter(rs.od + "/" + "fw.fastq"));
            BufferedWriter rwFastqWriter = new BufferedWriter(new FileWriter(rs.od + "/" + "rw.fastq"));

            while ((line = br.readLine()) != null) {
                if (header) {
                    header = false;
                } else {
                    entry = line.split("\t");   // 0: gene_id, 1: transcript_id, 2: count
                    Gene curGene = annotations.getGenes().get(entry[0]);
                    curStrand = curGene.getStrand();
                    Transcript curTranscript = annotations.getTranscriptById(entry[1]);
                    curTransSeq = rs.gu.getTransSeq(curGene, curTranscript, gse);
                    readcount = Integer.parseInt(entry[2]);

                    for (int i = 0; i < readcount; i++) {
                        fragmentData = rs.gu.getFragmentData(rs.frlength, rs.sd, rs.readlength, curTransSeq, fragmentData, curStrand);
                        if (!fragmentData.isInvalidSequence()) {
                            fragmentData.setCurGene(curGene);
                            fragmentData.settFwRegvec(rs.gu.getCoveredRegions(curTranscript, fragmentData, true, curTransSeq.length(), curStrand));
                            fragmentData.settRvRegvec(rs.gu.getCoveredRegions(curTranscript, fragmentData, false, curTransSeq.length(), curStrand));
                            mw.setFd(fragmentData);
                            mw.write(entry);
                            fwFastqWriter.write("@" + mw.getReadId() + "\n" + fragmentData.getFwRead() + "\n"
                                    + "+" + mw.getReadId() + "\n" + qualityScore + "\n");
                            rwFastqWriter.write("@" + mw.getReadId() + "\n" + fragmentData.getRvRead() + "\n"
                                    + "+" + mw.getReadId() + "\n" + qualityScore + "\n");
                            mw.incrementReadId();
                        } else {
                            fragmentData.setInvalidSequence(false);
                        }
                    }
                }
            }
            bwMappingInfo.close();
            fwFastqWriter.close();
            rwFastqWriter.close();
        }
    }
 */
}
