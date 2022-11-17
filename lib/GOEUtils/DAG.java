package GOEUtils;

import org.apache.commons.math3.distribution.HypergeometricDistribution;
import org.apache.commons.math3.stat.inference.KolmogorovSmirnovTest;

import java.util.*;

public class DAG {
    private HashMap<String, DAGNode> graph;
    private HashMap<String, EnrichmentEntry> enrichEntries;
    private HashMap<String, String> mappedGenes;
    private HashMap<String, EnrichmentEntry> annotatedEntries;
    private HashSet<String> trueEntries;
    private int deGenes;

    public DAG() {
        this.graph = new HashMap<>();
        this.enrichEntries = new HashMap<>();
        this.trueEntries = new HashSet<>();
        this.deGenes = 0;
        this.mappedGenes = new HashMap<>();
    }

    public HashMap<String, String> getMappedGenes() {
        return mappedGenes;
    }

    public void setMappedGenes(HashMap<String, String> mappedGenes) {
        this.mappedGenes = mappedGenes;
    }

    public HashSet<String> getTrueEntries() {
        return trueEntries;
    }

    public void setTrueEntries(HashSet<String> trueEntries) {
        this.trueEntries = trueEntries;
    }

    public int getDeGenes() {
        if (deGenes == 0) {
            for (Map.Entry<String, EnrichmentEntry> entry : enrichEntries.entrySet()) {
                EnrichmentEntry curEntry = entry.getValue();
                if (curEntry.isSignif() && curEntry.annotated) deGenes++;
            }
        }
        return deGenes;
    }

    public HashMap<String, EnrichmentEntry> getAnnotatedEntries() {
        if (annotatedEntries == null) {
            annotatedEntries = new HashMap<>();
            for (Map.Entry<String, EnrichmentEntry> entry : enrichEntries.entrySet()) {
                EnrichmentEntry curEntry = entry.getValue();
                if (curEntry.annotated) annotatedEntries.put(curEntry.getId(), curEntry);
            }
        }
        return annotatedEntries;
    }

    public void calcFdrs(ArrayList<DAGNode> nodes, String method) {
        if (method.equals("hg")) nodes.sort(Comparator.comparing(node -> node.getHgPval()));
        if (method.equals("fej")) nodes.sort(Comparator.comparing(node -> node.getFejPval()));
        if (method.equals("ks")) nodes.sort(Comparator.comparing(node -> node.getKsPval()));

        double[] allPvals = new double[nodes.size()];
        int y = 0;
        for (DAGNode node : nodes) {
            if (method.equals("hg")) allPvals[y] = node.getHgPval();
            if (method.equals("fej")) allPvals[y] = node.getFejPval();
            if (method.equals("ks")) allPvals[y] = node.getKsPval();
            y++;
        }
        double[] adjustedPvals = new double[allPvals.length];
        int m = allPvals.length;
        for (int i = m - 1; i >= 0; i--) {
            if (i == m - 1) {
                adjustedPvals[i] = allPvals[i];
            } else {
                double unadjustedPvalue = allPvals[i];
                int divideByM = i + 1;
                double left = adjustedPvals[i + 1];
                double right = (m / (double) divideByM) * unadjustedPvalue;
                adjustedPvals[i] = Math.min(left, right);
            }
        }
        int i = 0;
        for (DAGNode node : nodes) {
            if (method.equals("hg")) node.setHgFdr(adjustedPvals[i]);
            if (method.equals("fej")) node.setFejFdr(adjustedPvals[i]);
            if (method.equals("ks")) node.setKsFdr(adjustedPvals[i]);
            i++;
        }
    }

    public double hgPval(DAGNode node, int setSize, int overlap) {
        int totalGenes = getAnnotatedEntries().size();
        int deGenes = getDeGenes();
        HypergeometricDistribution hg = new HypergeometricDistribution(totalGenes, deGenes, setSize);
        return hg.upperCumulativeProbability(overlap);
    }

    public double fejPval(DAGNode node, int setSize, int overlap) {
        int totalGenes = getAnnotatedEntries().size() - 1;
        int deGenes = getDeGenes() - 1;
        HypergeometricDistribution hg = new HypergeometricDistribution(totalGenes, deGenes, setSize - 1);
        return hg.upperCumulativeProbability(overlap - 1);
    }

    public double[] ksMeasures(DAGNode node) {
        ArrayList<Double> setDistAl = new ArrayList<>();
        ArrayList<Double> brDistAl = new ArrayList<>();

        for (Map.Entry<String, EnrichmentEntry> gene : enrichEntries.entrySet()) {
            EnrichmentEntry curGene = gene.getValue();
            if (curGene.annotated) {
                if (!node.getGenes().contains(curGene.getId())) {
                    brDistAl.add(curGene.getFc());
                } else {
                    setDistAl.add(curGene.getFc());
                }
            }
        }
        KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();
        double[] setDist = new double[setDistAl.size()];
        for (int i = 0; i < setDistAl.size(); i++) {
            setDist[i] = setDistAl.get(i);
        }
        double[] brDist = new double[brDistAl.size()];
        for (int i = 0; i < brDistAl.size(); i++) {
            brDist[i] = brDistAl.get(i);
        }
        double ksPval = ks.kolmogorovSmirnovTest(setDist, brDist);
        double ksStat = ks.kolmogorovSmirnovStatistic(setDist, brDist);
        double[] measures = {ksPval, ksStat};
        return measures;
    }

    public int measuredAssocGenes(DAGNode node) {   // size column
        int size = 0;
        for (String assocGene : node.getGenes()) {
            if (enrichEntries.containsKey(assocGene)) {
                size++;
            }
        }
        return size;
    }

    // TODO combine upper method and this method

    public int getNoverlap(DAGNode node) {
        int noverlap = 0;
        for (String assocGene : node.getGenes()) {
            if (enrichEntries.containsKey(assocGene)) {
                EnrichmentEntry curGene = enrichEntries.get(assocGene);
                if (curGene.isSignif()) {
                    noverlap++;
                }
            }
        }
        return noverlap;
    }

    public String shortestPathToATrue(DAGNode node) {
        if (node.isTrue()) {
            return "";
        }
        int shortestPath = Integer.MAX_VALUE;
        ArrayList<String> finalPath = new ArrayList<>();

        for (String trueNode : trueEntries) {
            DAGNode curTrue = graph.get(trueNode);
            ArrayList<String> curPath = shortestPathBetween3(node, curTrue);
            if (curPath.size() < shortestPath) {
                shortestPath = curPath.size();
                finalPath = curPath;
            }
        }

        // build actual path
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < finalPath.size() - 1; i++) {
            sb.append(finalPath.get(i));
            sb.append("|");
        }
        sb.append(finalPath.get(finalPath.size() - 1)); // last entry
        return sb.toString();
    }

    public ArrayList<String> shortestPathBetween3(DAGNode cNode, DAGNode cTrue) {
        HashSet<String> commonAncestors = new HashSet<>(cNode.getAncestors());
        commonAncestors.retainAll(cTrue.getAncestors());
        commonAncestors.add(cNode.getId());
        ArrayList<String> singleNodePath = new ArrayList<>();
        singleNodePath.add(cNode.getId());
        cNode.getPaths().add(singleNodePath);

        commonAncestors.add(cTrue.getId());
        ArrayList<String> singleTruePath = new ArrayList<>();
        singleTruePath.add(cTrue.getId());
        cTrue.getPaths().add(singleTruePath);
        ArrayList<String> cShortest = new ArrayList<>();

        int cShortestSize = graph.size();    // better than Integer.MAXVALUE
        for (String ca : commonAncestors) {  // common ancestor
            ArrayList<String> nodeShortestPath = new ArrayList<>();
            int nodeSPathSize = graph.size();
            ArrayList<String> trueShortestPath = new ArrayList<>();
            int trueSPathSize = graph.size();

            if (ca.equals(cNode.getId())) {
                nodeShortestPath.add(cNode.getName() + " * ");
                nodeSPathSize = nodeShortestPath.size();
            }

            if (ca.equals(cTrue.getId())) {
                trueShortestPath.add(cTrue.getName() + " * ");
                trueSPathSize = trueShortestPath.size();
            }

            if (nodeShortestPath.size() != 1) {
                for (ArrayList<String> path : cNode.getPaths()) {
                    ArrayList<String> tempPath = new ArrayList<>();
                    tempPath.add(cNode.getName());
                    for (String term : path) {
                        DAGNode curTermNode = graph.get(term);
                        if (term.equals(ca)) {
                            tempPath.add(curTermNode.getName() + " * ");
                            if (tempPath.size() < nodeSPathSize) {
                                nodeShortestPath = tempPath;
                                nodeSPathSize = nodeShortestPath.size();
                                tempPath = new ArrayList<>();
                            } else {
                                break;
                            }
                        } else {
                            tempPath.add(curTermNode.getName()); // still going to ca
                        }
                    }
                }
            }

            if (trueShortestPath.size() != 1) {
                for (ArrayList<String> path : cTrue.getPaths()) {
                    ArrayList<String> tempPath = new ArrayList<>();
                    tempPath.add(cNode.getName());
                    for (String term : path) {
                        DAGNode curTermNode = graph.get(term);
                        if (term.equals(ca)) {
                            tempPath.add(curTermNode.getName() + " * ");
                            if (tempPath.size() < trueSPathSize) {
                                trueShortestPath = tempPath;
                                trueSPathSize = trueShortestPath.size();
                                tempPath = new ArrayList<>();
                            } else {
                                break;
                            }
                        } else {
                            tempPath.add(curTermNode.getName());
                        }
                    }
                }
            }
            if (nodeSPathSize + trueSPathSize < cShortestSize) {
                cShortestSize = nodeSPathSize + trueSPathSize;
                cShortest = nodeShortestPath;
                Collections.reverse(trueShortestPath);
                for (int i = 1; i < trueShortestPath.size(); i++) {
                    cShortest.add(trueShortestPath.get(i));
                }
            }
        }
        return cShortest;
    }

    public ArrayList<String> shortestPathBetween2(DAGNode node, DAGNode trueNode) {
        HashSet<String> commonAncestors = new HashSet<>(node.getAncestors());
        commonAncestors.retainAll(trueNode.getAncestors());

        ArrayList<String> finalPath = new ArrayList<>();
        int shortestPathSize = Integer.MAX_VALUE;

        // if directly connected
        if (node.getParents().contains(trueNode.getId()) ||
                trueNode.getParents().contains(node.getId())) {
            ArrayList<String> path = new ArrayList<>();
            path.add(node.getName());
            path.add(trueNode.getName());
            return path;

        } else {
            for (String ca : commonAncestors) {  // common ancestor
                ArrayList<String> nodePath = new ArrayList<>();
                int nodePathSize = graph.size();
                ArrayList<String> truePath = new ArrayList<>();
                int truePathSize = graph.size();

                for (ArrayList<String> path : node.getPaths()) {
                    ArrayList<String> potNodePath = new ArrayList<>();
                    potNodePath.add(node.getName());
                    for (String term : path) {
                        DAGNode curTerm = graph.get(term);
                        potNodePath.add(curTerm.getName());
                        if (term.equals(ca)) {
                            if (potNodePath.size() < nodePathSize) {
                                nodePath = potNodePath;
                                String lca = nodePath.get(nodePath.size() - 1);
                                lca += " * ";
                                nodePath.set(nodePath.size() - 1, lca);
                                nodePathSize = nodePath.size();
                            }
                        }
                    }
                }
                for (ArrayList<String> path : trueNode.getPaths()) {
                    ArrayList<String> potNodePath = new ArrayList<>();
                    potNodePath.add(trueNode.getName());
                    for (String term : path) {
                        DAGNode curTerm = graph.get(term);
                        potNodePath.add(curTerm.getName());
                        if (term.equals(ca)) {
                            if (potNodePath.size() < truePathSize) {
                                truePath = potNodePath;
                                truePathSize = truePath.size();
                            }
                        }
                    }
                }
                // TODO Arraylist collections.reverse(AL)
                if (nodePathSize + truePathSize - 1 < shortestPathSize) {   // subtract ca, since it is accounted for twice
                    shortestPathSize = nodePathSize + truePathSize - 1;
                    for (int i = truePathSize - 2; i >= 0; i--) {
                        nodePath.add(truePath.get(i));
                    }
                }
                finalPath = nodePath;
            }
            return finalPath;
        }
    }

    public int shortestPathBetween(DAGNode fop, DAGNode sop) {
        HashSet<String> commonAncestors = new HashSet<>(fop.getAncestors());
        commonAncestors.retainAll(sop.getAncestors());

        if (fop.getParents().contains(sop.getId())
                || sop.getParents().contains(fop.getId())) {
            return 1;
        } else {
            int shortestPath = graph.size();    // better than Integer.MAXVALUE
            for (String ca : commonAncestors) {  // common ancestor
                int fopMinDist = graph.size();
                int sopMinDist = graph.size();
                for (ArrayList<String> path : fop.getPaths()) {
                    for (String term : path) {
                        if (term.equals(ca)) {
                            if (path.indexOf(term) + 1 < fopMinDist) {
                                fopMinDist = path.indexOf(term) + 1;
                            }
                        } else if (term.equals(sop.getId())) {
                            if (path.indexOf(term) < fopMinDist) {
                                fopMinDist = path.indexOf(term);
                            }
                        }
                    }
                }
                for (ArrayList<String> path : sop.getPaths()) {
                    for (String term : path) {
                        if (term.equals(ca)) {
                            if (path.indexOf(term) + 1 < sopMinDist) {
                                sopMinDist = path.indexOf(term) + 1;
                            }
                        } else if (term.equals(fop.getId())) {
                            if (path.indexOf(term) < sopMinDist) {
                                sopMinDist = path.indexOf(term);
                            }
                        }
                    }
                }
                if (fopMinDist + sopMinDist < shortestPath) {
                    shortestPath = fopMinDist + sopMinDist;
                }
            }
            return shortestPath;
        }
    }

    public void fillNodePaths(DAGNode curNode, DAGNode curParent, ArrayList<String> prior) {
        for (String parent : curParent.getParents()) {
            var thisPath = new ArrayList<String>(prior);
            thisPath.add(parent);
            DAGNode parentNode = graph.get(parent);
            if (parentNode.getParents().size() == 0) {
                curNode.getPaths().add(thisPath);
            } else {
                fillNodePaths(curNode, parentNode, thisPath);
            }
        }
    }

    public void fillAllPaths() {
        for (Map.Entry<String, DAGNode> node : graph.entrySet()) {
            DAGNode curNode = node.getValue();
            fillNodePaths(curNode, curNode, new ArrayList<>());
        }
    }

    public void associatedGenes() {
        for (Map.Entry<String, DAGNode> node : graph.entrySet()) {
            DAGNode curNode = node.getValue();
            for (String ancestor : curNode.getAncestors()) {
                DAGNode curAncestor = graph.get(ancestor);
                curAncestor.getGenes().addAll(curNode.getGenes());
            }
        }
    }

    public void fillAncestors(DAGNode node, HashSet<String> nodeAncestors) {
        for (String parent : node.getParents()) {
            nodeAncestors.add(parent);
            DAGNode curParent = graph.get(parent);
            nodeAncestors.addAll(curParent.getParents());
            fillAncestors(curParent, nodeAncestors);
        }
    }

    public void propagateAncestors() {
        for (Map.Entry<String, DAGNode> node : graph.entrySet()) {
            DAGNode curNode = node.getValue();
            fillAncestors(curNode, curNode.getAncestors());
        }
    }


    public boolean isRelative(DAGNode fop, DAGNode sop) {
        if (fop.getAncestors().contains(sop.getId()) ||
                sop.getAncestors().contains(fop.getId())) {
            return true;
        }
        return false;
    }

    public int overlapping(DAGNode fop, DAGNode sop) {
        HashSet<String> overlap = new HashSet<>(fop.getGenes());
        overlap.retainAll(sop.getGenes());
        return overlap.size();
    }

    public double maxOverlap(DAGNode fop, DAGNode sop, int overlapsize) {
        double fopOverlap = (double) overlapsize / fop.getGenes().size();
        double sopOverlap = (double) overlapsize / sop.getGenes().size();
        double max = Math.max(fopOverlap, sopOverlap) * 100;
        return (Math.round(max * 100.0) / 100.0);
    }

    public DAGNode getNodeById(String id) {
        return graph.get(id);
    }

    public HashMap<String, DAGNode> getGraph() {
        return graph;
    }

    public HashMap<String, EnrichmentEntry> getEnrichEntries() {
        return enrichEntries;
    }

    public void setEnrichEntries(HashMap<String, EnrichmentEntry> enrichEntries) {
        this.enrichEntries = enrichEntries;
    }
}
