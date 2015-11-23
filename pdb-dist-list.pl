#!/usr/bin/perl -w
######################################################################################
#
# Ludo 24/02/2010
######################################################################################

use strict;

use constant USAGE =><<END;

SYNOPSIS:

 pdb-dist.pl pdb pair.dat

DESCRIPTION:

Prints the distance between atoms listed in pair.dat of pdb pdb.

END


if ($#ARGV!=1) {die USAGE;}

###################################################################################
my $pdb=$ARGV[0];
my $pair=$ARGV[1];
my $chainformat=0;		#if the pdb contain the chain put 1 otherwise leave 0
###################################################################################
# INIT
###################################################################################
my ($line,$nco,$nf,$i,$j,$c);
my ($d,$d1,$r1,$r2,$s,$rab,$frame,$atid);
my (%Dab,%co);
my %atm;
my @field;

##################################### DETERMINE IF CHAIN IS PRESENT ##############
$j=0;
open(f1,"$pdb")  or die "no pdb : $pdb";
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

################################################ READ PDB FILE ############################
open(f1,"$pdb")  or die "no pdb : $pdb";
while($line=<f1>){
    @field=split(" ",$line);
    if (index($field[0],"ATOM")==0) {
        $atid=$field[1];
        $atm{$atid}{atype}=$field[2];
        $atm{$atid}{restype}=$field[3];
        $atm{$atid}{resid}=$field[4];
        $atm{$atid}{x}=$field[5+$chainformat];
        $atm{$atid}{y}=$field[6+$chainformat];
        $atm{$atid}{z}=$field[7+$chainformat];
    }
}
close(f1);

################################################ READ PAIR LIST ############################
print "#list of contacts created with pdb-dist-list.pl\n#chainformat=$chainformat\n#pdb=$pdb\n#pair=$pair\n";
print "#ID   \tATID TYP RES:ID  \tATID TYP RES:ID  \tdist    \t|ri-rj|\n";
$i=0;
$j=0;
$nco=1;
open(f2,"$pair")  or die "no pair : $pair";
while($line=<f2>){
    @field=split(" ",$line);
    if (index($field[0],"#")!=0) {
        $i=$field[0];
        $j=$field[1];

        $d1=&Dist($atm{$i},$atm{$j});
        $r1=$atm{$i}{resid};
        $r2=$atm{$j}{resid};
        printf("%-5d\t%5d %3s %3s %-5d\t%4d %3s %3s %-5d\t%lf\t%d\n",$nco,$i,$atm{$i}{atype},$atm{$i}{restype},$r1,$j,$atm{$j}{atype},$atm{$j}{restype},$r2,$d1,Abs($r2-$r1));
        $nco++;
    }
}
close(f2);


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

