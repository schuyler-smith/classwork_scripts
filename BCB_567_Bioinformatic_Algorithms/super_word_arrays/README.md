# SWA.py #

This script is used to do a Super Word Alignment of a single nucleotide sequence.

## SWA.py
To optimally run the pipeline from the terminal, use the following command:  

```
python <SWA.py> <sequence.fasta> <word_model.txt> wlcut
```

## Input File(s)
Input sequence file should be in the fasta format 

```
>read_name
GATCGTAGAGTGAGACCTAGTGTTTG
```

Input word model file should be in the format 

	10111


as an example
### PARAMETER INPUT VALUES
The following are the parameters required by the script:

	wlcut	|	a positive integer

## Output

The ouput is printed to the standard output (ie to the screen of your terminal). 

To save the output to a file use the bash redirection operator to save to a file

```
pipeline_command > <output_file.txt>
```


## Author Information

Schuyler Smith
Ph.D. Student in Bioinformatics and Computational Biology
Department of Agriculural and Biosystems Engineering
Iowa State University  Ames, IA

e-mail @ schuyler.smith@iastate.edu




