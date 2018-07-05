import os
import sys
import subprocess
from Bio import SeqIO



subprocess.check_output(['Rscript', 'functions/thing.R'] + ['arg1', 'arg2', '...'])

