#!/bin/zsh

set -eu -o pipefail


PREFIX=main


compile_main () {
    clear
    xelatex $PREFIX.tex
}


compile_refs () {
    clear
    bibtex $PREFIX
}


# Compile the new files.
compile_main
compile_refs
compile_main
compile_main

