# SiteNuclearity

The **SiteNuclearity** is used to characterize metal sites in metalloprotein structures.
Specifically, it can be used to:

- Determine the first coordination sphere of metals, based on a distance
threshold.

- Determine which metals belong to the same site (i.e., those that share a
ligand or fall within a distance of 5 Å) and consequently identify
multinuclear sites, if present.

- Draw a graph representing the site and its first coordination sphere (as
shown in the figure).

<img width="582" height="453" alt="BEF_BE_130_B_3259" src="https://github.com/user-attachments/assets/6ba7c3db-7bc3-4f8d-b563-c28ba3362285" />

---

## Usage

To use **SiteNuclearity** you must have Java version 17 or greater.

Many Linux/Mac systems already have a version of Java installed. Check the installed version by running

```bash
 java -version
 ```

If the installed version is lower than version 17, you can still install version 17 as several versions can safely co-exist.
In that case, you can download the compressed tar file from [https://adoptium.net/temurin/releases/](https://adoptium.net/temurin/releases/) site and unpack it in any folder, e.g. /opt and run **SiteNuclearity** with

```bash
/opt/<java-version>/bin/java -jar SiteNuclearity-1.0.jar
```

Running **SiteNuclearity** without options will show a help with the various options available

```bash
usage: SiteNuclearity
 -d,--minimal_distance <arg>   Minimal distance for metals in the same
                               site. Default is 5.0 A.
 -e,--excluded_donors <arg>    Chemical symbols of the atoms (separated by
                               commas) excluded from metal ligands.
                               Default is C and H.
 -m,--metal <arg>              Chemical symbol of the metal of interest.
                               Default is all metals.
 -o,--overwrite                Overwrite existing files and directories.
 -p,--pdb <arg>                Local input PDB file.
 -t,--threshold <arg>          Coordination distance threshold. Default is
                               2.8 A.
 -w,--workdir <arg>            Directory where to find the input PDB files
                               and to write outputs. Default is ./
```

An example of use is

```bash
/opt/<java-version>/bin/java -jar SiteNuclearity-1.0.jar -p 12ca.pdb 
```
