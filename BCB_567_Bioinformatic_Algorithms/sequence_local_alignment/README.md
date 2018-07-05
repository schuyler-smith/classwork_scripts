# SDS_BLAST #

This script is used to do a basic local alignment (BLAST) of 2 nucleotide sequences.

## SDS_BLAST.sh
To optimally run the pipeline from the terminal, use the following command:  

```
sh SDS_BLAST.sh <sequence_1.fasta> <sequence_2.fasta> mismatch_score gap_open_penalty gap_extension_penalty"
```
The bash script is used to execute the pipeline scripts efficiently. It is possible to execute the pipeline without using the bash script and instead calling the python script (shown in the next section).

The python scripts used in the pipeline require [NumPy][NP]. If you do not have this package and plan to use python often I would recommend looking into [Anaconda][ANA], which includes this and many other useful python packages.

The pipeline as-is does not benefit from the bash wrapper. Many will argue that it is unneccessary in general, Python can do nearly anything that can be done in the bash. This is true but they are not done with the same efficiency. To allow my scripts scalability I alwyas use a bash wrapper, but also write them to be able to be called independantly.

To execute the pipeline without the use of the bash script, it takes a similar input to the command line:
```
python Scripts/blast.py <sequence_1.fasta> <sequence_2.fasta> mismatch_score gap_open_penalty gap_extension_penalty"
```

[NP]:https://www.scipy.org/scipylib/download.html
[ANA]:https://www.continuum.io/downloads
## Input File(s)
Input file(s) should be in the fasta format 
```
<read_name>
>GATCGTAGAGTGAGACCTAGTGTTTG
```
### PARAMETER INPUT VALUES
The following are the parameters required by the script:

		mismatch				|	a negative integer
		gap open penalty		|	non-negative integer
		gap extension penalty	|	positive integer

The base match score is always = 10 for this algorithm.

## Output

The ouput is printed to the standard output (ie to the screen of your terminal). To save the output to a file use the bash redirection operator to save to a file>

```
pipeline_command > <output_file.blast>
```


## Author Information

Schuyler Smith
Ph.D. Student in Bioinformatics and Computational Biology
Department of Agriculural and Biosystems Engineering
Iowa State University  Ames, IA

e-mail @ schuyler.smith@iastate.edu




