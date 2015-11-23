#!/bin/bash
###########################################################################################
#
# Script to generate 2 contact map CVs
# to use with PLUMED 2.1 given 2 initial PDBs
# 
# v1.1
#
# L. Sutto 1 Oct 2014
# l.sutto@ucl.ac.uk
###########################################################################################


###########################################################################################
###########################################################################################
#
#  ADAPT THE FOLLOWING PARAMETERS
#
###########################################################################################
###########################################################################################

# name of the 2 pdbs corresponding to the 2 different structures
# NB: they must have the same numbering (and number of atoms!)
pdb1="open-em.pdb"
pdb2="close-em.pdb"

# choose the identifiers name (e.g. active/inactive, open/close..)
name1="active"
name2="inactive"

# output directory (either relative to current or absolute)
outdir="out"

# consider contacts only involving atoms of this region
# defined with atom number (using the pdb numbering)
aastart=3287
aaend=3798


# below this line are usually good default values

# weight of regular contacts vs Salt Bridge contacts
weight_regular=1
weight_sb=3

# consider contacts only those pairs with:
# distance_pdb1 <= cutoff_in AND distance_pdb2 > cutoff_out
#                             OR
# distance_pdb2 <= cutoff_in AND distance_pdb1 > cutoff_out
# unit: Angstrom
cutoff_in=5
cutoff_out=7

# N and M exponent defining the smooth contact function (see PLUMED manual)
CMAP_N=6
CMAP_M=10

###########################################################################################
###########################################################################################
###########################################################################################


# USUALLY NO NEED TO TOUCH BELOW THIS LINE


###########################################################################################
# INITIALIZATION
. functions.sh
mkdir -p $outdir
checkexist $pdb1 $pdb2 pick-pairs.pl pdb-dist-list.pl salt-bridge.awk
###########################################################################################


out_step1="$outdir/co-$name1-$name2-CACBO-spec_$aastart-$aaend.list"
###########################################################################################
# CREATE specific contacts list (CA,CB,O atoms only in [aastart;aaend])
# within cutoff_in Ang in one conformation but more than cutoff_in Ang in the other
###########################################################################################

# pick the pairs unique to one or the other pdb
echo "1a. finding specific contacts pdb1.."
./pick-pairs.pl $pdb1 0 $cutoff_in 2 $pdb2 $cutoff_out 100.0 > $outdir/co-$name1-$cutoff_in-$cutoff_out.dat
checkerr
echo "1b. finding specific contacts pdb2.."
./pick-pairs.pl $pdb2 0 $cutoff_in 2 $pdb1 $cutoff_out 100.0 > $outdir/co-$name2-$cutoff_in-$cutoff_out.dat
checkerr

# filter out only those with an atom in [aastart;aaend] and of type CA, CB, O
grep -v "#" $outdir/co-$name1-$cutoff_in-$cutoff_out.dat  | awk -v s=$aastart -v e=$aaend 'function good(a){return ((a=="CA")||(a=="CB")||(a=="O"))} {if (good($3) && (($2>=s && $2<=e)||($6>=s && $6<=e)) && good($7) ) print $2,$6,$10}' > $outdir/co-$name1-CACBO-spec_$aastart-$aaend.list
checkerr
grep -v "#" $outdir/co-$name2-$cutoff_in-$cutoff_out.dat  | awk -v s=$aastart -v e=$aaend 'function good(a){return ((a=="CA")||(a=="CB")||(a=="O"))} {if (good($3) && (($2>=s && $2<=e)||($6>=s && $6<=e)) && good($7) ) print $2,$6,$10}' > $outdir/co-$name2-CACBO-spec_$aastart-$aaend.list
checkerr

# concatenate the 2 sets
# -> out_step1: ai aj dist_ref
nco1=`wc -l $outdir/co-$name1-CACBO-spec_$aastart-$aaend.list|awk '{print $1}'`
nco2=`wc -l $outdir/co-$name2-CACBO-spec_$aastart-$aaend.list|awk '{print $1}'`
echo "--> specific contacts of pdb1: $nco1"
echo "--> specific contacts of pdb2: $nco2"
cat $outdir/co-$name1-CACBO-spec_$aastart-$aaend.list $outdir/co-$name2-CACBO-spec_$aastart-$aaend.list > $out_step1
checkerr
echo "--> result: $out_step1"
echo


out_step2="$outdir/allsb.list"
###########################################################################################
# CREATE SaltBridges specific list: allsb.list
###########################################################################################

#generate a list of salt-bridge contacts (contacts between charged residues)
#nb set noh=0 to consider also H atoms.
echo "2. finding Salt Bridges.."
./pick-pairs.pl $pdb1 0 5 3 > $outdir/co-$name1.dat 
checkerr
./pick-pairs.pl $pdb2 0 5 3 > $outdir/co-$name2.dat 
checkerr

# useful for later to have also the distance!
./salt-bridge.awk VERBOSE=0 $outdir/co-$name1.dat > $outdir/sb-$name1-aiajd.list 
checkerr
./salt-bridge.awk VERBOSE=0 $outdir/co-$name2.dat > $outdir/sb-$name2-aiajd.list
checkerr

# concatenate them, filtering out all the duplicated SB (uses uniq -u command)
# -> out_step2: ai aj dist_ref
cat $outdir/sb-$name1-aiajd.list $outdir/sb-$name2-aiajd.list | sort -n -k1,2|awk '{print $3,$1,$2}' | uniq -u -f 1|awk '{print $2,$3,$1}' > $out_step2
checkerr

nco=`wc -l $out_step2|awk '{print $1}'`
echo "--> total Salt bridge contacts: $nco"
echo "--> result: $out_step2"
echo



###########################################################################################
# CREATE PLUMED2 CMAP INPUT FILES
###########################################################################################


# calculate the distance of the contact list of step1 for pdb1 and pdb2
# -> paste them so to have: ai aj d_pdb1 d_pdb2 d_ref
./pdb-dist-list.pl $pdb1 $out_step1 | grep -v "#" | awk '{print $2,$6,$10}' > $outdir/_tmp11
checkerr
./pdb-dist-list.pl $pdb2 $out_step1 | grep -v "#" | awk '{print $2,$6,$10}' > $outdir/_tmp12
checkerr
paste $outdir/_tmp11 $outdir/_tmp12 $out_step1 | awk '{print $1,$2,$3,$6,$9}' > $outdir/_tmp112
checkerr

# calculate the distance of the contact list of step2 for pdb1 and pdb2
# -> paste them so to have: ai aj d_pdb1 d_pdb2 d_ref
./pdb-dist-list.pl $pdb1 $out_step2 | grep -v "#" | awk '{print $2,$6,$10}' > $outdir/_tmp21
checkerr
./pdb-dist-list.pl $pdb2 $out_step2 | grep -v "#" | awk '{print $2,$6,$10}' > $outdir/_tmp22
checkerr
paste $outdir/_tmp21 $outdir/_tmp22 $out_step2 | awk '{print $1,$2,$3,$6,$9}'> $outdir/_tmp212
checkerr

# concatenate and use the weight:
# -> final-co.list: ai aj d_pdb1(nm) ref1 d_pdb2(nm) ref2 d_ref(nm) weight
# NB: the reference value must not be multiplied by the weight!
cat $outdir/_tmp112 | awk -v w=$weight_regular -v n=$CMAP_N -v m=$CMAP_M 'function ref(d,r0){if(d==r0) {return n/m;} else {return (1-(d/r0)**n)/(1-(d/r0)**m);}} {d1=$3/10.;d2=$4/10.;r0=$5/10.;print $1,$2,d1,ref(d1,r0),d2,ref(d2,r0),r0,w}' > $outdir/final-co.list
checkerr
# NB: the reference value must not be multiplied by the weight!
cat $outdir/_tmp212 | awk -v w=$weight_sb -v n=$CMAP_N -v m=$CMAP_M 'function ref(d,r0){if(d==r0) {return n/m;} else {return (1-(d/r0)**n)/(1-(d/r0)**m);}} {d1=$3/10.;d2=$4/10.;r0=$5/10.;print $1,$2,d1,ref(d1,r0),d2,ref(d2,r0),r0,w}' >> $outdir/final-co.list
checkerr

rm -f $outdir/_tmp* 

# finally, create the plumed2 inputs
./aux-plumed2-cmap.sh $outdir/final-co.list $name1 $name2 $CMAP_N $CMAP_M
checkerr

echo "Successfully terminated."
echo 

# create a plumed.dat to test the CVs:
echo "INCLUDE FILE=plum2-CMAP-$name1.dat" > plumed.dat
echo "INCLUDE FILE=plum2-CMAP-$name2.dat" >> plumed.dat
echo "PRINT ARG=cmap_$name1,cmap_$name2 FILE=COLVAR" >> plumed.dat

cat $pdb1 $pdb2 > pdb-2frames-$name1-$name2.pdb

echo "************************************************"
echo "TO TEST THE CVs ON THE 2 FRAMES RUN THE COMMAND:"
echo
echo "plumed driver --plumed plumed.dat --mf_pdb pdb-2frames-$name1-$name2.pdb"
echo
echo "and check the output in COLVAR file, it should look like:"
echo
echo "#! FIELDS time cmap_active cmap_inactive"
echo " 0.000000 0.000000 97.278820"
echo " 1.000000 97.278865 0.000000"
echo 
