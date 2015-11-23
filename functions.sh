###########################################################################
# COMMON FUNCTIONS
###########################################################################
checkerr(){
    if [ $? -ne 0 ]
    then
        echo "> Errore ! <"
        exit 1
    fi
}

checkexist()
{
    local file
    for file in $@
    do
        if [ ! -e $file ]; then
            echo "ERROR : $file does not exist, aborting"
            exit 1
        fi
    done
}

checkstring()
{
    if [ -n  "$1" ]; then
        echo
        echo "### STOPPING EXECUTION : ###"
        cat _tmp
        echo "############################"
        echo
        exit
    fi
}

exclmin()
{
    # function to return the smallest number in the range [1;max]
    # excluding all the n successive 'bad' numbers:
    #
    # number=exclmin max bad1 bad2 bad3 ... badn
    #
    # I use ( ) around the function so to make every variable local

    (
    ARGV=($@)
    max=$1
    nexcl=$(($#-1))

    for i in `seq 1 1 $max`
    do
        ok=1
        for n in `seq 1 1 $nexcl`
        do
            toexclude=${ARGV[$n]}
            if [ $i -eq $toexclude ]; then ok=0;fi
        done
        if [ $ok -eq 1 ]; then break;fi
    done
    echo $i
    )
}


