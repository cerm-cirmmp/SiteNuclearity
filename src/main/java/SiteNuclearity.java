
import org.apache.commons.cli.*;
import org.apache.commons.io.output.NullPrintStream;
import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.GroupType;
import org.graphper.api.*;
import org.graphper.api.attributes.*;
import org.graphper.draw.ExecuteException;

import java.io.File;
import java.io.IOException;
import java.util.*;

import static java.lang.System.exit;

public class SiteNuclearity {

    private static String pdbFile = null;

    private static Double threshold = 2.8;
    private static Double minimal_distance = 5.0;
    private static String workdir = "." + File.separator;
    private static String metal = null;
    private static List<String> notDonors = Arrays.asList("H", "C");
    private static final Constants constants = new Constants();
    private static List<String> metals = constants.getMetals();

    private static void parseArgumentLine(String[] args) {

        Options opt = new Options();
        opt.addOption("p", "pdb", true, "Local input PDB file.");
        opt.addOption("e", "excluded_donors", true, "Chemical symbols of the atoms (separated by commas) excluded from metal ligands. Default is C and H.");
        opt.addOption("t", "threshold", true, "Coordination distance threshold. Default is 2.8 A.");
        opt.addOption("d", "minimal_distance", true, "Minimal distance for metals in the same site. Default is 5.0 A.");
        opt.addOption("o", "overwrite", false, "Overwrite existing files and directories.");
        opt.addOption("m", "metal", true, "Chemical symbol of the metal of interest. Default is all metals.");
        opt.addOption("w", "workdir", true, "Directory where to find the input PDB files and to write outputs. Default is ./");

        CommandLineParser parser = new DefaultParser();

        HelpFormatter formatter = new HelpFormatter();

        try {
            CommandLine cmd = parser.parse(opt, args);

            if (!cmd.hasOption("p")) {
                throw new ParseException("A PDB file is required.");
            }
            if (cmd.hasOption("p")) {
                pdbFile = cmd.getOptionValue("p");
            }

            if (cmd.hasOption("t")) {
                String tshold = cmd.getOptionValue("t");
                try {
                    threshold = Double.parseDouble(tshold);
                } catch (Exception e) {
                    throw new ParseException("Invalid threshold. This must be a number.");
                }
            }
            if (cmd.hasOption("d")) {
                String minDist = cmd.getOptionValue("d");
                try {
                    minimal_distance = Double.parseDouble(minDist);
                } catch (Exception e) {
                    throw new ParseException("Invalid distance. This must be a number.");
                }
            }

            if (cmd.hasOption("w")) {
                workdir = cmd.getOptionValue("w");
                File path = new File(workdir);
                /*if (fileExists(workdir) && !cmd.hasOption("o")) {
                    throw new ParseException(workdir + " is an existing directory... Remove/rename it or use the -o option");
                } else if (fileExists(workdir) && cmd.hasOption("o")) {
                    deleteDirectory(path);
                }
                createDir(path);*/
                if (!fileExists(workdir)) {
                    throw new ParseException(workdir + " doesn't exists...Set an existing directory where I can find pdb files.");
                }
            }

            if (cmd.hasOption("m")) {

                metal = cmd.getOptionValue("m").toUpperCase();
                if (!metals.contains(metal)) {
                    throw new ParseException(metal + " invalid metal.");
                }
                metals = List.of(metal);
            }

            if (cmd.hasOption("e")) {
                try {
                    notDonors = Arrays.asList(cmd.getOptionValue("e").replace(" ", "").split(","));
                } catch (Exception e) {
                    throw new ParseException("!!! Unexpected exception: " + e);
                }
            }

        } catch (ParseException exp) {
            printLine();
            System.out.println("!!! Unexpected exception: " + exp.getMessage());
            printLine();
            formatter.printHelp("SiteNuclearity", opt);
            exit(-1);
        }
    }

    public static void createDir(File path) {
        path.mkdir();
    }

    /**
     * Check if a file or directory exists
     *
     * @param filename the path to check
     * @return true if file/directory exists or false if not
     */
    public static Boolean fileExists(String filename) {
        File f = new File(filename);
        if (!f.exists()) {
            return false;
        }
        return true;
    }

    /**
     * remove a directory
     *
     * @param path the directory to remove
     * @return true if directory is removed
     */
    public static boolean deleteDirectory(File path) {
        if (path.exists()) {
            File[] files = path.listFiles();
            for (int i = 0; i < files.length; i++) {
                if (files[i].isDirectory()) {
                    deleteDirectory(files[i]);
                } else {
                    files[i].delete();
                }
            }
        }
        return (path.delete());
    }

    private static void printLine() {
        for (int i = 0; i < 81; i++) {
            System.out.print("-");
        }
        System.out.println();
    }


    // Crea un grafico per ogni metallo.
    // Ogni grafico è costruito legando il metallo ad ogni donatore e ogni donatore al suo residuo di appartenenza.
    // In questo metodo se due donatori appartengono allo stesso residuo sono comunque legati a due nodi diversi con lo stesso nome.
    private static void createGraph(List<List<Atom>> sites, LinkedHashMap<String, List<Atom>> donors, String workdir) {
        HashMap<String, String> donors_colors = new HashMap<String, String>() {
            {
                put("B", "#ffb5b5");
                put("C", "#33ff33");
                put("S", "#FFFF00");
                put("O", "#FF0000");
                put("P", "#ff7f00");
                put("I", "#940094");
                put("BR", "#a62929");
                put("CL", "#1ff01f");
                put("F", "#b3ffff");
                put("X", "#FFFFFF");
                put("SE", "#ffa100");

            }
        };

        for (int i = 0; i < sites.size(); i++) {
            File path = new File(workdir + "/site" + (i + 1) + "_" + sites.get(i).size());
            createDir(path);
            for (Atom a : sites.get(i)) {
                Graphviz.GraphvizBuilder builder = Graphviz.digraph();
                Node node1 = Node.builder().label(a.getElement().toString()).fillColor(Color.GREY).build();

                String key = a.getElement() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChainId() + "_" + a.getPDBserial();
                String grphFile = a.getElement() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChain().getName() + "_" + a.getPDBserial();
                List<Atom> res = donors.get(key);

                for (Atom atm : res) {
                    Node n1 = Node.builder().label(atm.getName()).fillColor(Color.BLUE).fontColor(Color.WHITE).build();
                    Node n2 = null;
                    String k = a.getElement() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChainId() + "_" + a.getPDBserial();
                    for (Atom donRes : donors.get(k)) {
                        if (donRes == atm) {
                            n2 = Node.builder().label(donRes.getGroup().getPDBName() + "_" +
                                    donRes.getGroup().getResidueNumber() + "(" + donRes.getGroup().getChain().getName() + ")").fillColor(Color.ORANGE).build();
                            break;
                        }
                    }

                    builder.addLine(node1, n1, n2);
                }
                try {
                    builder.build().toFile(FileType.PNG).save(path.toString(), grphFile);
                } catch (IOException | ExecuteException e) {
                    System.out.println(e);
                }
            }
        }

    }

    // Crea un grafico per ogni sito.
    // Ogni grafico è costruito legando i metalli ad ogni donatore e ogni donatore al suo residuo di appartenenza.
    // In questo metodo se due metalli o due donatori sono legati allo stesso donatore o appartengono allo stesso residuo, sono legati ad uno stesso nodo.
    private static void createGraph3(List<List<Atom>> sites, LinkedHashMap<String, List<Atom>> donors, String workdir) {
        for (int i = 0; i < sites.size(); i++) {
            File path = new File(workdir + "/site" + (i + 1) + "_" + sites.get(i).size());
            createDir(path);
            Graphviz.GraphvizBuilder builder = Graphviz.digraph();
            List<NodesDTO> nodesDTOS = new ArrayList<>();
            Map<String, Node> thirdRankNodeMap = new HashMap<>();
            Map<Atom, Node> secondRankNodeMap = new HashMap<>();

            for (Atom a : sites.get(i)) {

                NodesDTO nodesDTO = new NodesDTO();
                nodesDTO.setNode1(Node.builder().label(a.getElement().toString()).fillColor(Color.GREY).build());
                String resMetal = a.getGroup().getPDBName() + "_" + a.getGroup().getResidueNumber().toString() + "_" + a.getGroup().getChain().getName();
                String key = a.getElement() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChainId() + "_" + a.getPDBserial();
                List<Atom> res = donors.get(key);
                List<Node> level1 = new ArrayList<>();
                List<Node> level2 = new ArrayList<>();

                for (Atom atm : res) {
                    Node n1 = secondRankNodeMap.computeIfAbsent(atm, l -> Node.builder().label(atm.getName()).fillColor(Color.BLUE).fontColor(Color.WHITE).build());
                    level1.add(n1);

                    String k = a.getElement() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChainId() + "_" + a.getPDBserial();
                    for (Atom donRes : donors.get(k)) {
                        if (donRes == atm) {
                            Color color;
                            if(donRes.getGroup().getType() == GroupType.HETATM) {
                                String donMetal = donRes.getGroup().getPDBName() + "_" + donRes.getGroup().getResidueNumber().toString() + "_" + donRes.getGroup().getChain().getName();
                                if(resMetal.equals(donMetal)) {
                                    color = Color.PINK;
                                } else {
                                    color = Color.GOLD;
                                }
                            } else {
                                color = Color.ORANGE;
                            }
                            String label = donRes.getGroup().getPDBName() + "_" +
                                    donRes.getGroup().getResidueNumber() + "(" + donRes.getGroup().getChain().getName() + ")";
                            Node n2 = thirdRankNodeMap.computeIfAbsent(label, l -> Node.builder().label(label).fillColor(color).build());
                            level2.add(n2);
                            break;
                        }
                    }

                }
                nodesDTO.setLevel1(level1);
                nodesDTO.setLevel2(level2);
                nodesDTOS.add(nodesDTO);
            }
            for (int x = 0; x < nodesDTOS.size(); x++) {
                for (int y = 0; y < nodesDTOS.get(x).getLevel1().size(); y++) {
                    Node n1 = nodesDTOS.get(x).getNode1();
                    Node n2 = nodesDTOS.get(x).getLevel1().get(y);
                    Node n3 = nodesDTOS.get(x).getLevel2().get(y);
                    // Only one edge gonna be rendered if id is duplicated
                    Line line = Line.builder(n2, n3).id(n2.nodeAttrs().getLabel() + "^" + n3.nodeAttrs().getLabel()).build();
                    builder.addLine(n1, n2).addLine(line);
                }

            }
            try {
                builder.build().toFile(FileType.PNG).save(path.toString(), "" + (i+1));
            } catch (IOException | ExecuteException e) {
                System.out.println(e);
            }

        }
    }

    // Crea un grafico per ogni metallo.
    // Ogni grafico è costruito legando il metallo ad ogni donatore e ogni donatore al suo residuo di appartenenza.
    // In questo metodo se due donatori appartengono allo stesso residuo sono legati allo stesso nodo.
    private static void createGraph4(List<Atom> sites, LinkedHashMap<String, List<Atom>> donors, String workdir, int i) {

            File path = new File(workdir + "/site" + (i + 1) + "_" + sites.size());
            createDir(path);

            for (Atom a : sites) {
                Graphviz.GraphvizBuilder builder = Graphviz.digraph();
                String resMetal = a.getGroup().getPDBName() + "_" + a.getGroup().getResidueNumber().toString() + "_" + a.getGroup().getChain().getName();
                Map<String, Node> thirdRankNodeMap = new HashMap<>();
                Map<Atom, Node> secondRankNodeMap = new HashMap<>();
                List<NodesDTO> nodesDTOS = new ArrayList<>();
                NodesDTO nodesDTO = new NodesDTO();
                nodesDTO.setNode1(Node.builder().label(a.getElement().toString()).fillColor(Color.GREY).build());

                String key = a.getElement() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChainId() + "_" + a.getPDBserial();

                String grphFile = a.getGroup().getPDBName() + "_" + a.getElement().toString().toUpperCase() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChain().getName() + "_" + a.getPDBserial();
                List<Atom> res = donors.get(key);
                List<Node> level1 = new ArrayList<>();
                List<Node> level2 = new ArrayList<>();

                for (Atom atm : res) {
                    Node n1 = secondRankNodeMap.computeIfAbsent(atm, l -> Node.builder().label(atm.getName()).fillColor(Color.BLUE).fontColor(Color.WHITE).build());
                    level1.add(n1);

                    String k = a.getElement() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChainId() + "_" + a.getPDBserial();
                    for (Atom donRes : donors.get(k)) {
                        if (donRes == atm) {
                            Color color;
                            if(donRes.getGroup().getType() == GroupType.HETATM) {
                                String donMetal = donRes.getGroup().getPDBName() + "_" + donRes.getGroup().getResidueNumber().toString() + "_" + donRes.getGroup().getChain().getName();
                                if(resMetal.equals(donMetal)) {
                                    color = Color.PINK;
                                } else {
                                    color = Color.GOLD;
                                }
                            } else {
                                color = Color.ORANGE;
                            }
                            String label = donRes.getGroup().getPDBName() + "_" +
                                    donRes.getGroup().getResidueNumber() + "(" + donRes.getGroup().getChain().getName() + ")";
                            Node n2 = thirdRankNodeMap.computeIfAbsent(label, l -> Node.builder().label(label).fillColor(color).build());
                            level2.add(n2);
                            break;
                        }
                    }

                }
                nodesDTO.setLevel1(level1);
                nodesDTO.setLevel2(level2);
                nodesDTOS.add(nodesDTO);

                for (int x = 0; x < nodesDTOS.size(); x++) {
                    for (int y = 0; y < nodesDTOS.get(x).getLevel1().size(); y++) {
                        Node n1 = nodesDTOS.get(x).getNode1();
                        Node n2 = nodesDTOS.get(x).getLevel1().get(y);
                        Node n3 = nodesDTOS.get(x).getLevel2().get(y);
                        // Only one edge gonna be rendered if id is duplicated
                        Line line = Line.builder(n2, n3).id(n2.nodeAttrs().getLabel() + "^" + n3.nodeAttrs().getLabel()).build();
                        builder.addLine(n1, n2).addLine(line);
                    }

                }
                try {
                    builder.build().toFile(FileType.PNG).save(path.toString(), grphFile);
                } catch (IOException | ExecuteException e) {
                    System.out.println(e);
                }
            }



    }

    // Create one graph for each site without bridged ligands
    private static void createGraph2(List<List<Atom>> sites, LinkedHashMap<String, List<Atom>> donors, String workdir) {

        for (int i = 0; i < sites.size(); i++) {

            File path = new File(workdir + "/site" + (i + 1) + "_" + sites.get(i).size());
            createDir(path);
            Graphviz.GraphvizBuilder builder = Graphviz.digraph();
            List<NodesDTO> nodesDTOS = new ArrayList<>();
            for (Atom a : sites.get(i)) {

                NodesDTO nodesDTO = new NodesDTO();
                nodesDTO.setNode1(Node.builder().label(a.getElement().toString()).fillColor(Color.GREY).build());

                String key = a.getElement() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChainId() + "_" + a.getPDBserial();
                List<Atom> res = donors.get(key);
                List<Node> level1 = new ArrayList<>();
                List<Node> level2 = new ArrayList<>();

                for (Atom atm : res) {
                    level1.add(Node.builder().label(atm.getName()).fillColor(Color.BLUE).fontColor(Color.WHITE).build());

                    String k = a.getElement() + "_" + a.getGroup().getResidueNumber() + "_" + a.getGroup().getChainId() + "_" + a.getPDBserial();
                    for (Atom donRes : donors.get(k)) {
                        if (donRes == atm) {
                            level2.add(Node.builder().label(donRes.getGroup().getPDBName() + "_" +
                                    donRes.getGroup().getResidueNumber() + "(" + donRes.getGroup().getChain().getName() + ")").fillColor(Color.ORANGE).build());
                            break;
                        }
                    }

                }
                nodesDTO.setLevel1(level1);
                nodesDTO.setLevel2(level2);
                nodesDTOS.add(nodesDTO);
            }
            for (int x = 0; x < nodesDTOS.size(); x++) {
                for (int y = 0; y < nodesDTOS.get(x).getLevel1().size(); y++) {
                    builder.addLine(nodesDTOS.get(x).getNode1(), nodesDTOS.get(x).getLevel1().get(y), nodesDTOS.get(x).getLevel2().get(y));
                }

            }
            try {
                builder.build().toFile(FileType.PNG).save(path.toString(), "" + (i+1));
            } catch (IOException | ExecuteException e) {
                System.out.println(e);
            }

        }

    }


    public static LinkedHashMap<String, List<Atom>> copy(
            LinkedHashMap<String, List<Atom>> original) {
        LinkedHashMap<String, List<Atom>> copy = new LinkedHashMap<>();
        for (java.util.Map.Entry<String, List<Atom>> entry : original.entrySet()) {
            copy.put(entry.getKey(),
                    new ArrayList<>(entry.getValue()));
        }
        return copy;
    }


    public static void main(String[] args) {
        System.setErr(NullPrintStream.INSTANCE);

        parseArgumentLine(args);
        printLine();
        System.out.println("Input parameters");
        printLine();
        System.out.println("PDB: " + pdbFile);
        System.out.println("Threshold: " + threshold);
        System.out.println("Minimal Distance: " + minimal_distance);
        System.out.println("Metal: " + ((metal == null) ? "ALL" : metal));
        System.out.println("Workdir: " + workdir);
        System.out.println("Not Donors: " + notDonors.toString());
        printLine();

        PDBMethods pdbMethods = new PDBMethods();
        ReturnDTO metalsInPDB = pdbMethods.findMetal(pdbFile, workdir, metals);

        // System.out.println("Metals: " + metalsInPDB.getMetalsInPDB());

        LinkedHashMap<String, List<Atom>> donors = pdbMethods.findDonors(metalsInPDB.getMetalsInPDB(),
                metalsInPDB.getStructure().getChains(), threshold, notDonors);


        LinkedHashMap<String, List<Atom>> donorsFinale = copy(donors);

        //System.out.println("Donors :" + donors);
        //System.out.println("DonorsFinale :" + donorsFinale);

        List<List<Atom>> sites = pdbMethods.findSites2(metalsInPDB.getMetalsInPDB(), donors, minimal_distance);

        //System.out.println("Sites :" + sites);

        //createGraph(sites, donorsFinale, workdir);
        for (int i = 0; i < sites.size(); i++) {
                createGraph4(sites.get(i), donorsFinale, workdir, i);
        }
        createGraph3(sites, donorsFinale, workdir);
    }
}
