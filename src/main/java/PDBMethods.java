import org.biojava.nbio.structure.*;
import org.biojava.nbio.structure.align.util.AtomCache;
import org.biojava.nbio.structure.io.PDBFileReader;
import org.biojava.nbio.structure.io.StructureFiletype;

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

public class PDBMethods {

    public ReturnDTO findMetal(String pdbFilename, String pdbPath, List<String> cofactor) {
        List<String> metals;
        if (!cofactor.isEmpty()) {
            metals = cofactor;
        } else {
            Constants constants = new Constants();
            metals = constants.getMetals();
        }
        ReturnDTO returnDTO = new ReturnDTO();

        AtomCache cache = new AtomCache();
        cache.setPath(pdbPath + File.pathSeparator);
        cache.setFiletype(StructureFiletype.CIF);
        Path path = Paths.get(pdbPath, pdbFilename);
        StructureIO.setAtomCache(cache);
        Structure structure = null;
        try {
            structure = StructureIO.getStructure(path.toString());

        } catch (Exception e) {
            System.out.println(e);
            System.out.println("Unable to read PDB file!");
            return returnDTO;
        }

        List<Atom> metalsInPDB = new ArrayList<Atom>();
        List<Atom> metalsAltConf = new ArrayList<>();
        for (Chain c : structure.getChains()) {
            for (Group g : c.getAtomGroups(GroupType.HETATM)) {
                for (Atom a : g.getAtoms()) {
                    if (metals.contains(a.getElement().toString().toUpperCase())) {
                        //System.out.println(a);
                        metalsInPDB.add(a);

                    }
                }
                if (g.hasAltLoc()) {
                    for (Group altLoc: g.getAltLocs()) {
                        for (Atom t : altLoc.getAtoms()) {
                            if (metals.contains(t.getElement().toString().toUpperCase())) {
                                metalsAltConf.add(t);
                            }
                        }
                    }
                }
            }
        }

        returnDTO.setMetalsAlcConf(metalsAltConf);
        returnDTO.setMetalsInPDB(metalsInPDB);
        returnDTO.setStructure(structure);
        return returnDTO;
    }

    public LinkedHashMap<String, List<Atom>> findDonors(List<Atom> metals, List<Chain> chains, Double threshold, List<String> notDonors) {

        LinkedHashMap<String, List<Atom>> metal2donors = new LinkedHashMap<>();
        for (Atom m : metals) {
            List<String> resNum = new ArrayList<>();
            String key = m.getElement() + "_" + m.getGroup().getResidueNumber() + "_" + m.getGroup().getChainId() + "_" + m.getPDBserial();

            List<Atom> donors = new ArrayList<>();

            for (Chain c : chains) {
                //System.out.println(c.getAtomSequence());
                for (Group g : c.getAtomGroups()) {

                    //System.out.println(g.getChain().getAtomSequence());
                    for (Atom a : g.getAtoms()) {

                        if (g.getType() == GroupType.HETATM && metals.contains(a)) {
                            continue;
                        }
                        double x = m.getX() - a.getX();
                        double y = m.getY() - a.getY();
                        double z = m.getZ() - a.getZ();
                        //System.out.println(a.getGroup().getChain().getAtomSequence());
                        double distance = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2));
                        if (distance < threshold && !a.getName().equals(m.getName()) && !notDonors.contains(a.getElement().toString().toUpperCase()) ) {
                            resNum.add(g.getResidueNumber() + "_" + g.getChain());
                            donors.add(a);
                        }
                    }
                    metal2donors.put(key, donors);
                }
            }

            if (metal2donors.get(key).isEmpty()) {
                System.out.printf("No ligands found for metal %s (try to change distance threshold)", m.getName());
            }
        }

        //System.out.println("metal2donors: " + metal2donors);
        return metal2donors;
    }

    public List<List<Atom>> findSites2(List<Atom> metals, LinkedHashMap<String, List<Atom>> donors, Double minDistance) {
        //System.out.println("DonorsPrima: " + donors);
        List<List<Atom>> metalSites = new ArrayList<>();

        LinkedHashMap<String, List<String>> donorsTemp = new LinkedHashMap<>();

        for (Map.Entry<String, List<Atom>> set : donors.entrySet()) {
            List<String> dtemp = new ArrayList<>();

            for (Atom a: set.getValue()) {
                dtemp.add(a.getGroup().getPDBName() + "_" + a.getGroup().getChainId() + "_" + a.getGroup().getResidueNumber());
            }
            donorsTemp.put(set.getKey(), dtemp);
        }

        //System.out.println("donorsTemp: " + donorsTemp);
        //System.out.println("metals: " + metals);
        while (!metals.isEmpty()) {
            List<Atom> metalsAssigned = new ArrayList<>();
            List<Atom> tempList = new ArrayList<>();
            String keyTemp = metals.getFirst().getElement() + "_" + metals.getFirst().getGroup().getResidueNumber() + "_" + metals.getFirst().getGroup().getChainId() + "_" + metals.getFirst().getPDBserial();

            tempList.add(metals.getFirst());
            List<Atom> ligTemp = donors.get(keyTemp);
            metalsAssigned.add(metals.getFirst());
            metals.removeFirst();

            boolean new_member = true;
            while (new_member) {
                new_member = false;
                for (int i = 0; i < metals.size(); i++) {
                    String key1 = metals.get(i).getElement() + "_" + metals.get(i).getGroup().getResidueNumber() + "_" + metals.get(i).getGroup().getChainId() + "_" + metals.get(i).getPDBserial();

                    for (Atom mt : metalsAssigned) {
                        String key2 = mt.getElement() + "_" + mt.getGroup().getResidueNumber() + "_" + mt.getGroup().getChainId() + "_" + mt.getPDBserial();

                        Set<String> commons = donorsTemp.get(key1).stream()
                                .distinct()
                                .filter(donorsTemp.get(key2)::contains)
                                .collect(Collectors.toSet());

                        double x = metals.get(i).getX() - mt.getX();
                        double y = metals.get(i).getY() - mt.getY();
                        double z = metals.get(i).getZ() - mt.getZ();
                        //System.out.println(a.getGroup().getChain().getAtomSequence());
                        double distance = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2));
                        if (!commons.isEmpty() || distance < minDistance) {
                            tempList.add(metals.get(i));
                            metalsAssigned.add(metals.get(i));
                            metals.remove(i);

                            // Attenzione! Questo ciclo modifica anche la LinkedHashMap donors
                            // aggiungendo di fatto dei ligands alla lista già in essere
                            for (Atom ligands : donors.get(key1)) {
                                if (!ligTemp.contains(ligands)) {
                                    ligTemp.add(ligands);
                                }
                            }

                            new_member = true;
                            break;
                        }
                    }

                }
            }

            metalSites.add(tempList);

        }
        //System.out.println("DonorsTemp: " + donorsTemp);
        //System.out.println("Donors: " + donors);
        return metalSites;
    }

    public LinkedHashMap<String, List<List<String>>> findSites(List<Atom> metals, LinkedHashMap<String, List<List<String>>> donors, Double minDistance) {
        LinkedHashMap<String, List<List<String>>> metalSites = new LinkedHashMap<>();

        LinkedHashMap<String, List<String>> donorsTemp = new LinkedHashMap<>();
        for (Map.Entry<String, List<List<String>>> set : donors.entrySet()) {
            List<String> dtemp = new ArrayList<>();

            for (int i = 0; i < set.getValue().size(); i++) {
                dtemp.add(set.getValue().get(i).get(1) + "_" + set.getValue().get(i).get(2) + "_" + set.getValue().get(i).get(4));
            }
            donorsTemp.put(set.getKey(), dtemp);
        }

        //System.out.println("donorsTemp: " + donorsTemp);
        while (!metals.isEmpty()) {
            List<Atom> metalsAssigned = new ArrayList<>();
            String keyTemp = metals.getFirst().getName() + "_" + metals.getFirst().getGroup().getResidueNumber() + "_" + metals.getFirst().getGroup().getChain().getName();
            List<List<String>> ligTemp = donors.get(keyTemp);
            metalsAssigned.add(metals.getFirst());
            metals.removeFirst();

            boolean new_member = true;
            while (new_member) {
                new_member = false;
                for (int i = 0; i < metals.size(); i++) {
                    String key1 = metals.get(i).getName() + "_" + metals.get(i).getGroup().getResidueNumber() + "_" + metals.get(i).getGroup().getChain().getName();

                    for (Atom mt : metalsAssigned) {
                        String key2 = mt.getName() + "_" + mt.getGroup().getResidueNumber() + "_" + mt.getGroup().getChain().getName();

                        Set<String> commons = donorsTemp.get(key1).stream()
                                .distinct()
                                .filter(donorsTemp.get(key2)::contains)
                                .collect(Collectors.toSet());
                        //System.out.println("commons: " + commons);
                        double x = metals.get(i).getX() - mt.getX();
                        double y = metals.get(i).getY() - mt.getY();
                        double z = metals.get(i).getZ() - mt.getZ();
                        //System.out.println(a.getGroup().getChain().getAtomSequence());
                        double distance = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2) + Math.pow(z, 2));
                        if (!commons.isEmpty() || distance < minDistance) {
                            keyTemp += ";" + key1;
                            metalsAssigned.add(metals.get(i));
                            metals.remove(i);

                            // Attenzione! Questo ciclo modifica anche la LinkedHashMap donors
                            // aggiungendo di fatto dei ligands alla lista già in essere
                            for (List<String> ligands : donors.get(key1)) {
                                if (!ligTemp.contains(ligands)) {
                                    ligTemp.add(ligands);
                                }
                            }

                            new_member = true;
                            break;
                        }
                    }

                }
            }
            metalSites.put(keyTemp, ligTemp);

        }
        // System.out.println("metalSites: " + metalSites);
        return metalSites;
    }
}
