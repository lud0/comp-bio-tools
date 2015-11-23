#!/usr/bin/awk -f 

function charged(res){
    if (res=="ARG" || res=="LYS") return 1;
    if (res=="ASP" || res=="GLU") return -1; 
    return 0;
}
function okatom(aa,res) {
    if (res=="ARG" && aa=="CZ") return 1;
    if (res=="LYS" && aa=="NZ") return 1;
    if (res=="ASP" && aa=="CG") return 1;
    if (res=="GLU" && aa=="CD") return 1;
    return 0;
}
BEGIN {
    VERBOSE=1

    if (ARGC==1) {
        print "usage: salt-bridge.awk [VERBOSE=0/1]  nco.dat";
        print
        print "SB definition: contact between"
        print " CZ_ARG(+) -- CG_ASP(-)"
        print " CZ_ARG(+) -- CD_GLU(-)"
        print " NZ_LYS(+) -- CG_ASP(-)"
        print " NZ_LYS(+) -- CD_GLU(-)"
        print
        print "Prints the salt bridges found in the contact list file nco.dat"
        print "Expecting the column format:"
        print
        print "ID        AT1  ATM RES RESID  AT2 ATM RES RESID  dist"
        print "8771      4449 OC2 GLN 274    1765 OE1 GLN 110   4.654825"
        print
        print "if variable VERBOSE=0 ony prints [AT1 AT2 dist] of the salt bridge otherwise more info"

        exit 0;
    }
}
{
    if (VERBOSE==1) {
        if (index($1,"#")!=1 && charged($4)*charged($8)==-1 && okatom($3,$4) && okatom($7,$8)) print $2,$3,$4,$5,$6,$7,$8,$9
    } else {
        if (index($1,"#")!=1 && charged($4)*charged($8)==-1 && okatom($3,$4) && okatom($7,$8)) print $2,$6,$10
    }
}
