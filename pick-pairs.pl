#!/usr/bin/perl -w
######################################################################################
#
# Ludo 25/10/2013
######################################################################################

use strict;

use constant USAGE =><<END;

SYNOPSIS:

 pick-pairs.pl pdb1 min1 max1 bbd [pdb2 min2 max2]

DESCRIPTION:

Prints the atoms whose distance in pdb1 is in the range [min1;max1] belonging to residues i,j such that j >= i+bbd.
If pdb2 is present, then the same pair has also to be in the range [min2;max2] in pdb2.
Excluding hydrogen atoms if \$noh is set 1.


END


if ($#ARGV!=3 && $#ARGV!=6) {die USAGE;}

###################################################################################
my $pdb1=$ARGV[0];
my $min=$ARGV[1];
my $max=$ARGV[2];
my $bbd=$ARGV[3];
my ($pdb2,$min2,$max2,$duppdb);
if ($#ARGV!=6) {
    $pdb2=$ARGV[4];
    $min2=$ARGV[5];
    $max2=$ARGV[6];
    $duppdb=1;
} else {
    $pdb2="";
    $min2=-999999;
    $max2=999999;
    $duppdb=0;
}
my $noh=1;			    # set to 0 to include hydrogen
my $debug=0;

###################################################################################
# INIT
###################################################################################
my ($line,$nco,$nf,$i,$j,$c,$chainformat);
my ($d,$d2,$s,$rab,$frame,$atype);
my (%Dab,%co);
my (%atm,%atm2);
my @field;

##################################### DETERMINE IF CHAIN IS PRESENT ##############
$j=0;
open(f1,"$pdb1")  or die "no pdb : $pdb1";
while($line=<f1>){
    @field=split(" ",$line);
    if (index($field[0],"ATOM")==0) {
        if ($field[4] =~ /[A-Za-z]/) {
            # field[4] is alphabetic
            $chainformat=1;
        } else {
            $chainformat=0;
        }
        last;
    }
}
close(f1);
####################################### READ PDB FILE ############################
$i=0;
$j=0;
open(f1,"$pdb1")  or die "no pdb : $pdb1";
while($line=<f1>){
    @field=split(" ",$line);
    if (index($field[0],"ATOM")==0) {
        $atm{$i}{$j}{atid}=$field[1];
        $atm{$i}{$j}{atype}=$field[2];
        $atm{$i}{$j}{restype}=$field[3];
        $atm{$i}{$j}{resid}=$field[4+$chainformat];
        $atm{$i}{$j}{x}=$field[5+$chainformat];
        $atm{$i}{$j}{y}=$field[6+$chainformat];
        $atm{$i}{$j}{z}=$field[7+$chainformat];
        $atype=$field[2];

        #chose only heavy atoms if noh
        if (($noh) && (index($atype,"H")==0 || $atype =~ /^\dH/ )) {
            if ($debug) {print "#excluding $field[1] $atype\n";};
            $j--;
        }
        $j++;
    }
}
$atm{$i}{natom}=$j;
close(f1);


if ($duppdb==1) {
    $i=0;
    $j=0;
    open(f1,"$pdb2")  or die "no pdb : $pdb2";
    while($line=<f1>){
        @field=split(" ",$line);
        if (index($field[0],"ATOM")==0) {
            $atm2{$i}{$j}{x}=$field[5+$chainformat];
            $atm2{$i}{$j}{y}=$field[6+$chainformat];
            $atm2{$i}{$j}{z}=$field[7+$chainformat];
            $atype=$field[2];

            #chose only heavy atoms if noh
            if (($noh) && (index($atype,"H")==0 || $atype =~ /^\dH/ )) {
                if ($debug) {print "#excluding $field[1] $atype\n";};
                $j--;
            }
            $j++;
        }
    }
    close(f1);
}

print "#list of contacts created with pick-pairs.pl\n#chainformat=$chainformat\n#noh=$noh\n#min=$min\n#max=$max\n#min2=$min2\n#max2=$max2\n#bbd=$bbd\n#pdb=$pdb1\n#pdb2=$pdb2\n#totatom=$j\n";
print "#ID   \tATID TYP RES:ID  \tATID TYP RES:ID  \tdist    \tdist2   \t|ri-rj|\n";
my $d1;
my $r1;
my $r2;

$nco=0;
for ($i=0;$i<$atm{0}{natom};$i++) {
    for ($j=$i+1;$j<$atm{0}{natom};$j++) {
        $d1=&Dist($atm{0}{$i},$atm{0}{$j});
        if ($duppdb==1) {
            # in case of double PDB, enforce also the second condition
            $d2=&Dist($atm2{0}{$i},$atm2{0}{$j});
        } else {
            $d2=0;
        }
        $r1=$atm{0}{$i}{resid};
        $r2=$atm{0}{$j}{resid};
        if (($d1>$min)&&($d1<$max)&&($r2>=$r1+$bbd)&&($d2>$min2)&&($d2<$max2) ) { 
            $nco++;
            printf("%-5d\t%5d %3s %3s %-5d\t%4d %3s %3s %-5d\t%lf\t%lf\t%d\n",$nco,$atm{0}{$i}{atid},$atm{0}{$i}{atype},$atm{0}{$i}{restype},$atm{0}{$i}{resid},$atm{0}{$j}{atid},$atm{0}{$j}{atype},$atm{0}{$j}{restype},$atm{0}{$j}{resid},$d1,$d2,Abs($r2-$r1));
        }
    }
}



###################################################################################
# FUNCTIONS
###################################################################################
sub Dist
{
    my $atm1 = shift;
    my $atm2 = shift;

    my $d;

    $d=sqrt(($atm1->{x} - $atm2->{x})**2 + ($atm1->{y} - $atm2->{y})**2 + ($atm1->{z} - $atm2->{z})**2);
    return $d;
}

sub Abs
{
    my $i = shift;

    if ($i<0) {return -$i;}
    return $i;
}

