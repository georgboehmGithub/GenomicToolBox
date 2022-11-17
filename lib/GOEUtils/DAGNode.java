package GOEUtils;

import java.util.ArrayList;
import java.util.HashSet;

public class DAGNode {
    private final String id;
    private final HashSet<String> parents;
    private boolean isTrue;
    private final String name;
    private HashSet<String> ancestors;
    private HashSet<String> genes;
    private ArrayList<ArrayList<String>> paths;
    private String shortestToTrue;
    private int noverlap;
    private int size;
    private double hgPval;
    private double hgFdr;
    private double fejPval;
    private double fejFdr;
    private double ksPval;
    private double ksFdr;
    private double ksStat;

    public DAGNode(String id, HashSet<String> parents, String name) {
        this.id = id;
        this.name = name;
        this.parents = parents;
        this.genes = new HashSet<>();
        this.ancestors = new HashSet<>();
        this.paths = new ArrayList<>();
        this.isTrue = false;
        this.shortestToTrue = "";
    }

    public double getHgFdr() {
        return hgFdr;
    }

    public void setHgFdr(double hgFdr) {
        this.hgFdr = hgFdr;
    }

    public double getFejFdr() {
        return fejFdr;
    }

    public void setFejFdr(double fejFdr) {
        this.fejFdr = fejFdr;
    }

    public double getKsFdr() {
        return ksFdr;
    }

    public void setKsFdr(double ksFdr) {
        this.ksFdr = ksFdr;
    }

    public double getHgPval() {
        return hgPval;
    }

    public void setHgPval(double hgPval) {
        this.hgPval = hgPval;
    }

    public double getFejPval() {
        return fejPval;
    }

    public void setFejPval(double fejPval) {
        this.fejPval = fejPval;
    }

    public double getKsPval() {
        return ksPval;
    }

    public void setKsPval(double ksPval) {
        this.ksPval = ksPval;
    }

    public double getKsStat() {
        return ksStat;
    }

    public void setKsStat(double ksStat) {
        this.ksStat = ksStat;
    }

    public int getNoverlap() {
        return noverlap;
    }

    public void setNoverlap(int noverlap) {
        this.noverlap = noverlap;
    }

    public int getSize() {
        return size;
    }

    public void setSize(int size) {
        this.size = size;
    }

    public String getName() {
        return name;
    }

    public String getShortestToTrue() {
        return shortestToTrue;
    }

    public void setShortestToTrue(String shortestToTrue) {
        this.shortestToTrue = shortestToTrue;
    }

    public ArrayList<ArrayList<String>> getPaths() {
        return paths;
    }

    public void setPaths(ArrayList<ArrayList<String>> paths) {
        this.paths = paths;
    }

    public String getId() {
        return id;
    }

    public HashSet<String> getParents() {
        return parents;
    }

    public HashSet<String> getGenes() {
        if (genes == null) {
            this.genes = new HashSet<>();
        }
        return genes;
    }

    public boolean isTrue() {
        return isTrue;
    }

    public void setTrue(boolean aTrue) {
        isTrue = aTrue;
    }

    public void setGenes(HashSet<String> genes) {
        this.genes = genes;
    }

    public HashSet<String> getAncestors() {
        return ancestors;
    }

    public void setAncestors(HashSet<String> ancestors) {
        this.ancestors = ancestors;
    }
}
