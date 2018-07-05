#!/bin/bash

the_path=$(dirname "${0}")
SEQ_1=${1}
SEQ_2=${2}
MISMATCH=${3}
GAP_OPEN=${4}
GAP_EXTEND=${5}

python $the_path/Scripts/blast.py $SEQ_1 $SEQ_2 $MISMATCH $GAP_OPEN $GAP_EXTEND