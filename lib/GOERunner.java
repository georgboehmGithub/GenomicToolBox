import GOEUtils.DAG;
import GOEUtils.DAGNode;
import GOEUtils.EnrichmentEntry;
import com.beust.jcommander.JCommander;

import java.io.*;
import java.util.*;
import java.util.zip.GZIPInputStream;

public class GOERunner {
    private File obo;
    private File mapping;
    private int minsize;
    private int maxsize;
    private String root;
    private String enrich;
    private String overlapout;
    private String mappingtype;
    private String outputPath;
    private DAG gos;

    private static GOERunner runner = new GOERunner();

    private GOERunner() {
    }

    public static GOERunner getInstance() {
        return runner;
    }

    protected static void initRunner(Args argsParser) {
        runner.obo = new File(argsParser.obo);
        runner.mapping = new File(argsParser.mapping);
        runner.minsize = Integer.parseInt(argsParser.minsize);
        runner.maxsize = Integer.parseInt(argsParser.maxsize);
        runner.root = argsParser.root;
        runner.enrich = argsParser.enrich;
        runner.overlapout = argsParser.overlapout;
        runner.mappingtype = argsParser.mappingtype;
        runner.outputPath = argsParser.outputPath;
        runner.gos = new DAG();
    }

    protected static void parseGo() throws IOException {
        String id;
        String modifier;
        String assocGo;
        HashMap<String, DAGNode> graph = runner.gos.getGraph();
        try (BufferedReader br = new BufferedReader(new InputStreamReader
                (new GZIPInputStream(new FileInputStream(runner.mapping))))) {
            String line;
            String[] lineSeparated;
            while ((line = br.readLine()) != null) {
                if (line.charAt(0) != '!') {
                    lineSeparated = line.split("\t");
                    id = lineSeparated[2];
                    modifier = lineSeparated[3];
                    assocGo = lineSeparated[4];
                    if (modifier.equals("")) {  // no association qualifier modifier.
                        if (graph.containsKey(assocGo)) {
                            runner.gos.getMappedGenes().put(id, id);
                            graph.get(assocGo).getGenes().add(id);      //TODO something is wrong here
                        }
                    }
                }
            }
        }
    }

    protected static void parseEnsembl() {
        String id = "";
        HashMap<String, DAGNode> graph = runner.gos.getGraph();
        boolean first = false;
        try (BufferedReader br = new BufferedReader(new FileReader(runner.mapping))) {
            String line;
            String[] lineSeparated;
            while ((line = br.readLine()) != null) {
                if (!first) {   // Skip header
                    first = true;
                } else {
                    lineSeparated = line.split("\t");
                    id = lineSeparated[1];
                    if (id.equals("")) {
                        continue;
                    }
                    String[] allGos = lineSeparated[2].split("\\|");
                    for (String go : allGos) {
                        if (graph.containsKey(go)) {
                            runner.gos.getMappedGenes().put(id, id);
                            graph.get(go).getGenes().add(id);
                        }
                    }
                }
            }
        } catch (IOException e) {
            System.err.println("Error parsing ensembl file. Make sure mappingtype = 'ensembl' and file" +
                    "format is correct.");
        }
    }

    protected static void parseObo() {
        boolean newTerm = false;
        String id = "";
        String[] nameParts;
        StringBuilder name = new StringBuilder();
        var parents = new HashSet<String>();
        HashMap<String, DAGNode> graph = runner.gos.getGraph();
        try (BufferedReader br = new BufferedReader(new FileReader(runner.obo))) {
            String line;
            while ((line = br.readLine()) != null) {
                if (!newTerm) {
                    if (line.startsWith("[")) {
                        newTerm = true;
                    }
                } else {
                    if (line.equals("")) { // term end
                        DAGNode newGO = new DAGNode(id, new HashSet<>(parents), name.toString());
                        graph.put(id, newGO);
                        newTerm = false;
                        parents.clear();
                        name = new StringBuilder();
                    } else if (line.startsWith("name:")) {
                        nameParts = line.split(" "); // contains parts of name separated by whitespace
                        for (int i = 1; i < nameParts.length; i++) {
                            name.append(nameParts[i]);
                            if (i != nameParts.length - 1) {
                                name.append(" ");
                            }
                        }
                    } else if (line.startsWith("is_obsolete")) {
                        newTerm = false;    // ignore this term
                        name = new StringBuilder();
                    } else if (line.startsWith("id")) {
                        id = line.split(" ")[1];
                    } else if (line.startsWith("namespace")) {
                        if (!line.split(" ")[1].equals(runner.root)) {
                            newTerm = false;
                            name = new StringBuilder();
                        }
                    } else if (line.startsWith("is_a")) {
                        parents.add(line.split(" ")[1]);
                    }
                }
            }
        } catch (IOException e) {
            System.err.println("Error parsing obo file");
        }
    }

    protected static void parseMappingFile() throws IOException {
        boolean validInput = false;
        if (runner.mappingtype.equals("go")) {
            validInput = true;
            parseGo();
        } else if (runner.mappingtype.equals("ensembl")) {
            validInput = true;
            parseEnsembl();
        } else {
            System.err.println("Mappingtype should be either 'go' or 'ensembl'");
        }
    }

    protected static void parseEnrichFile() {
        try (BufferedReader br = new BufferedReader(new FileReader(runner.enrich))) {
            String line;
            String[] lineSeparated;
            boolean first = false;
            while ((line = br.readLine()) != null) {
                if (!line.startsWith("#")) {    // non-simulated DAGNode
                    if (!first) {
                        first = true;
                    } else {
                        lineSeparated = line.split("\t");
                        String id = lineSeparated[0];
                        String fc = lineSeparated[1];
                        String signif = lineSeparated[2];
                        EnrichmentEntry newEntry = new EnrichmentEntry(id, fc, signif);
                        if (runner.gos.getMappedGenes().containsKey(id)) {
                            newEntry.annotated = true;
                        }
                        runner.gos.getEnrichEntries().put(id, newEntry);
                    }
                } else {    // simulated DAG Node
                    String goTerm = line.substring(1, line.length());
                    DAGNode simulTrue = runner.gos.getGraph().get(goTerm);
                    simulTrue.setTrue(true);
                    runner.gos.getTrueEntries().add(goTerm);
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) throws IOException {
        Args argsParser = new Args();
        JCommander jct = JCommander.newBuilder()
                .addObject(argsParser)
                .build();
        jct.parse(args);
        if (argsParser.isHelp()) {
            jct.usage();
        } else {
            initRunner(argsParser);
            parseObo();
            parseMappingFile();
            runner.gos.propagateAncestors();    // set ancestors via parents
            runner.gos.associatedGenes();   // give every ancestor of nodes their gene set
            runner.gos.fillAllPaths();  // set every path to root for each node

            parseEnrichFile();

            if (!runner.overlapout.equals("")) {    // only optional param
                BufferedWriter bw = new BufferedWriter(new FileWriter(runner.overlapout));
                bw.write("term1" + "\t" + "term2" + "\t" + "is_relative" + "\t" +
                        "path_length" + "\t" + "num_overlapping" + "\t" + "max_ov_percent" + "\n");
                DAGNode fop;    // first of pair
                DAGNode sop; // second of pair
                HashMap<String, Integer> solvedPairs = new HashMap<>();

                for (Map.Entry<String, DAGNode> fopTerm : runner.gos.getGraph().entrySet()) {
                    fop = fopTerm.getValue();
                    if (fop.getGenes().size() >= runner.minsize && fop.getGenes().size() <= runner.maxsize) {
                        for (Map.Entry<String, DAGNode> sopTerm : runner.gos.getGraph().entrySet()) {
                            sop = sopTerm.getValue();
                            if (!fop.getId().equals(sop.getId())) {
                                if (sop.getGenes().size() >= runner.minsize && sop.getGenes().size() <= runner.maxsize) {
                                    String fopKey = fop.getId() + "_" + sop.getId();
                                    String sopKey = sop.getId() + "_" + fop.getId();
                                    if (!solvedPairs.containsKey(fopKey) && !solvedPairs.containsKey(sopKey)) {
                                        solvedPairs.put(fopKey, 1);
                                        // calculations, otherwise pair skipped
                                        boolean isRelative = runner.gos.isRelative(fop, sop);
                                        int numOverlapping = runner.gos.overlapping(fop, sop);
                                        double maxOverlap = runner.gos.maxOverlap(fop, sop, numOverlapping);
                                        int shortestPath = runner.gos.shortestPathBetween(fop, sop);
                                        if (numOverlapping > 0) {
                                            bw.write(fop.getId() + "\t" + sop.getId() + "\t" + isRelative + "\t" + shortestPath
                                                    + "\t" + numOverlapping + "\t" + maxOverlap + "\n");
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                bw.close();
            }

            BufferedWriter enrichBw = new BufferedWriter(new FileWriter(runner.outputPath));
            enrichBw.write("term" + "\t" + "name" + "\t" + "size" + "\t" +
                    "is_true" + "\t" + "noverlap" + "\t" + "hg_pval" + "\t" +
                    "hg_fdr" + "\t" + "fej_pval" + "\t" + "fej_fdr" + "\t" + "ks_stat" + "\t" +
                    "ks_pval" + "\t" + "ks_fdr" + "\t" + "shortest_path_to_a_true" + "\n");

            DAGNode curTerm;
            for (Map.Entry<String, DAGNode> node : runner.gos.getGraph().entrySet()) {
                curTerm = node.getValue();
                if (curTerm.getGenes().size() >= runner.minsize && curTerm.getGenes().size() <= runner.maxsize) {
                    int size = runner.gos.measuredAssocGenes(curTerm);
                    int noverlap = runner.gos.getNoverlap(curTerm);
                    double hgPval = runner.gos.hgPval(curTerm, size, noverlap);
                    double fejPval = runner.gos.fejPval(curTerm, size, noverlap);
                    double[] ksMeasures = runner.gos.ksMeasures(curTerm);
                    double ksPval = ksMeasures[0];
                    double ksStat = ksMeasures[1];
                    String sptat = runner.gos.shortestPathToATrue(curTerm);
                    curTerm.setShortestToTrue(sptat);
                    curTerm.setSize(size);
                    curTerm.setNoverlap(noverlap);
                    curTerm.setHgPval(hgPval);
                    curTerm.setFejPval(fejPval);
                    curTerm.setKsPval(ksPval);
                    curTerm.setKsStat(ksStat);
                }
            }
            ArrayList<DAGNode> nodes = new ArrayList<>();
            for (Map.Entry<String, DAGNode> node : runner.gos.getGraph().entrySet()) {
                curTerm = node.getValue();
                if (curTerm.getGenes().size() >= runner.minsize && curTerm.getGenes().size() <= runner.maxsize) {
                    nodes.add(curTerm);
                }
            }
            runner.gos.calcFdrs(nodes, "hg");
            runner.gos.calcFdrs(nodes, "fej");
            runner.gos.calcFdrs(nodes, "ks");

            for (DAGNode node : nodes) {
                enrichBw.write(node.getId() + "\t" + node.getName() + "\t" + node.getSize() + "\t"
                        + node.isTrue() + "\t" + node.getNoverlap() + "\t" + node.getHgPval() + "\t"
                        + node.getHgFdr() + "\t" + node.getFejPval() + "\t" + node.getFejFdr() + "\t"
                        + node.getKsStat() + "\t" + node.getKsPval() + "\t" + node.getKsFdr() + "\t"
                        + node.getShortestToTrue() + "\n");
            }
            enrichBw.close();
        }
    }
}
