#!/bin/zsh

set -eu -o pipefail


clean_prefix ()
{
    # Unlink the output and temporary files starting with the prefix.
    prefix=$1
    rm -fv $prefix.aux $prefix.bbl $prefix.blg $prefix.lof $prefix.log $prefix.lot $prefix.nlo $prefix.out $prefix.pdf $prefix.synctex.gz $prefix.toc
}


clean_pdfa ()
{
    # Unlink files for PDF/A format.
    rm -fv $1.xmpdata creationdate.lua creationdate.timestamp pdfa.xmpi
}


# Remove all existing output and temporary files, if any.
for name in main abstract intro results strat discuss methods suppfigs
do
	echo "Cleaning files for $name"
	clean_prefix "$name"
done

echo "Cleaning files for PDF/A"
clean_pdfa main

