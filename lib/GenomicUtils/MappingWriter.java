package GenomicUtils;

import java.io.BufferedWriter;
import java.io.IOException;

public class MappingWriter {
    private FragmentData fd;
    private BufferedWriter bw;
    private int readId;

    public MappingWriter(FragmentData fd, BufferedWriter bw) {
        this.fd = fd;
        this.readId = 0;
        this.bw = bw;
    }

    public void write(String [] entry) throws IOException {
        this.bw.write(readId + "\t" + fd.getCurGene().getChromosome() + "\t" + entry[0] + "\t" + entry[1] + "\t" +
                fd.getfLocalRegion().getStart() + "-" + (fd.getfLocalRegion().getStart() + fd.getReadLength()) + "\t" +
                (fd.getfLocalRegion().getEnd() - fd.getReadLength()) + "-" + fd.getfLocalRegion().getEnd() + "\t" +
                fd.gettFwRegvec().toString() + "\t" + fd.gettRvRegvec().toString() + "\t" + fd.getFwMut() + "\t" +
                fd.getRvMut() + "\n");
    }

    public void incrementReadId() {
        this.readId += 1;
    }

    public int getReadId() {
        return this.readId;
    }

    public void setFd(FragmentData fd) {
        this.fd = fd;
    }
}
