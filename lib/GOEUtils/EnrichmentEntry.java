package GOEUtils;

public class EnrichmentEntry {
    private String id;
    private double fc;
    private boolean signif;
    public boolean annotated = false;

    public EnrichmentEntry(String id, String fc, String signif) {
        this.id = id;
        this.fc = Double.parseDouble(fc);
        this.signif = Boolean.parseBoolean(signif);
    }

    public double getFc() {
        return fc;
    }

    public void setFc(double fc) {
        this.fc = fc;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public boolean isSignif() {
        return signif;
    }

    public void setSignif(boolean signif) {
        this.signif = signif;
    }
}
