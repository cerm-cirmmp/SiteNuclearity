import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;

import java.util.List;

public class ReturnDTO {

    private List<Atom> metalsInPDB;

    private List<Atom> metalsAlcConf;

    Structure structure;

    public List<Atom> getMetalsInPDB() {
        return metalsInPDB;
    }

    public void setMetalsInPDB(List<Atom> metalsInPDB) {
        this.metalsInPDB = metalsInPDB;
    }

    public List<Atom> getMetalsAlcConf() {
        return metalsAlcConf;
    }

    public void setMetalsAlcConf(List<Atom> metalsAlcConf) {
        this.metalsAlcConf = metalsAlcConf;
    }

    public Structure getStructure() {
        return structure;
    }

    public void setStructure(Structure structure) {
        this.structure = structure;
    }
}
