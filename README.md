# comp-bio-tools
A repository of scripts for data analysis and manipulation useful for Molecular Dynamics simulations of proteins.

## Main scripts

* reweight.py: Standalone python script version of the Time-independent Free Energy reconstruction script (a.k.a. reweight) based on the algorithm proposed by [Tiwary and Parrinello JPCB 2014](http://pubs.acs.org/doi/abs/10.1021/jp504920s). The script is meant to be used as an analysis tool for a Molecular Dynamics simulation where the Metadynamics enhanced sampling technique is used to calculate a system's Free Energy.

* create-plumed2-cmaps.sh: Given 2 PDB files, creates the [PLUMED2](https://github.com/plumed/plumed2) input files using 2 contact maps (CMAPs) collective variables (CVs). The CMAPs are generated taking the atoms that are within a certain range in one pdb and outside another range in the second pdb. The typical situation is to use Metadynamics to enhance the sampling of a conformational transition between 2 known protein conformers.
This script depends on a number of other scripts in this same directory.


## Utilities - Miscellaneous

* aux-plumed2-cmap.sh: auxiliary script that creates 2 contact map files
* functions.sh: repository of commonly used BASH functions
* pdb-dist-list.pl: prints the distance between atoms listed in input file of given PDB file
* pick-pairs.pl: extract atom pairs within certain distance range in one or more PDB files
* salt-bridges.awk: given a list of atom pairs, returns the pairs that may be forming a salt bridge according to their residue type and atom type
