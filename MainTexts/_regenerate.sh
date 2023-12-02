#!/bin/zsh

set -eu -o pipefail


OUTPUT="maintext.pdf"

# Remove any old intermediate files.
rm -fv $OUTPUT
./_clean.sh

# Compile the latest files into a PDF.
./_compile.sh

# Rename the PDF from "main" to "thesis".
mv main.pdf $OUTPUT

# Open the new PDF.
open -a Preview $OUTPUT

# Remove the latest intermediate files.
./_clean.sh

