#!/bin/bash

# source common functions
. functions.sh

# check arguments
if [ $# -ne 5 ]
then
        echo "usage: aux-plumed2-cmap.sh MAPFILE name1 name2 N M"
        echo
        echo "format of MAPFILE: ai aj d_pdb1(nm) ref1 d_pdb2(nm) ref2 d_ref(nm) weight"
        exit
fi

map=$1
name1=$2
name2=$3
CMAP_N=$4
CMAP_M=$5

checkexist $map

out="plum2-CMAP-$name1.dat"

echo "CONTACTMAP ..." > $out
awk '{n+=1;printf("ATOMS%d=%d,%d REFERENCE%d=%f WEIGHT%d=%.2f\n",n,$1,$2,n,$4,n,$8)}' $map >> $out
awk -v nn=$CMAP_N -v mm=$CMAP_M '{n+=1;printf("SWITCH%d={RATIONAL R_0=%f D_0=0.0 NN=%d MM=%d}\n",n,$7,nn,mm)}' $map >> $out
echo "LABEL=cmap_$name1" >> $out
echo "CMDIST" >> $out
echo "... CONTACTMAP" >> $out

echo "Successfully generated the map:" $out

out="plum2-CMAP-$name2.dat"

echo "CONTACTMAP ..." > $out
awk '{n+=1;printf("ATOMS%d=%d,%d REFERENCE%d=%f WEIGHT%d=%.2f\n",n,$1,$2,n,$6,n,$8)}' $map >> $out
awk -v nn=$CMAP_N -v mm=$CMAP_M '{n+=1;printf("SWITCH%d={RATIONAL R_0=%f D_0=0.0 NN=%d MM=%d}\n",n,$7,nn,mm)}' $map >> $out
echo "LABEL=cmap_$name2" >> $out
echo "CMDIST" >> $out
echo "... CONTACTMAP" >> $out

echo "Successfully generated the map:" $out

