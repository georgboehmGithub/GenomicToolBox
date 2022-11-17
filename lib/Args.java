import com.beust.jcommander.Parameter;

import java.util.ArrayList;
import java.util.List;

public class Args {
    @Parameter(names = {"-maxsize"}, description = "maximal number of associated genes for GO entry")
    public String maxsize = "";

    @Parameter(names = {"-minsize"}, description = "minimal number of associated genes to GO entry")
    public String minsize = "";

    @Parameter(names = {"-enrich"}, description = "input file for the enrichment analysis")
    public String enrich = "";

    @Parameter(names = {"-overlapout"}, description = "output path for DAG entries shared mapped genes information")
    public String overlapout = "";

    @Parameter(names = {"-mappingtype"}, description = "format of how 'mapping' is provided")
    public String mappingtype = "";

    @Parameter(names = {"-mapping"}, description = "file which holds gene to GO mapping information")
    public String mapping = "";

    @Parameter(names = {"-root"}, description = "defines GO namespace")
    public String root = "";

    @Parameter(names = {"-obo"}, description = "input obo filepath")
    public String obo = "";

    @Parameter(names = {"-o"}, description = "Path to output directory")
    public String outputPath = "";

    @Parameter
    public List<String> parameters = new ArrayList<>();

    @Parameter(names = {"--help", "-h"}, help = true, description = "Show this message")
    private boolean help = false;

    public boolean isHelp() {
        return help;
    }
}
