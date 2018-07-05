# MultipleSequenceAlignment.py #

This script is used to do a Super Word Alignment of multiple nucleotide sequences.

## MSA.py
To optimally run the pipeline from the terminal, use the following command:  

```
python <MSA.py> <sequence.fasta> <word_model.txt> wlcut
```

## Input File(s)
Input sequence file should be in the fasta format 

```
>read_name_1
GATCGTAGAGTGAGACCTAGTGTTTG
>read_name_2
GATCGTAGAGTGAGACCTAGTGTTTG
```

Input word model file should be in the binary format 

	10111


as an example
### PARAMETER INPUT VALUES
The following are the parameters required by the script:

	wlcut	|	a positive integer

## Output

The ouput is printed to the standard output (ie to the screen of your terminal). 
```
>Seq1
ANACNCGNGNNTNTTNTN
>Seq2
NANACNCGNGNNNTNTTNT
>Seq3
NANACNCGNGNNTNTTNTNN

Word model: 101
wlcut: 2
The length of a longest chain of superword blocks: 15
The number of superword blocks in the chain: 3

Block 1:
Seq1       1        ANACNC
Seq2       2        ANACNC
Seq3       2        ANACNC

Block 2:
Seq1       4        CNCGNG
Seq2       5        CNCGNG
Seq3       5        CNCGNG

Block 3:
Seq1       12       TNTTNT
Seq2       14       TNTTNT
Seq3       13       TNTTNT
```

To save the output to a file use the bash redirection operator to save to a file

```
MSA_pipeline_command > <output_file.txt>
```


## Author Information

Schuyler Smith
Ph.D. Student in Bioinformatics and Computational Biology
Department of Agriculural and Biosystems Engineering
Iowa State University  Ames, IA

e-mail @ schuyler.smith@iastate.edu




