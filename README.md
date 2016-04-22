MAAPs
------

### Usage
```bash
maaps --fasta examples/Den4_MAAPS_TestData16.fasta
```

Maaps currently only takes one argument, `--fasta`. All other arguments are ignored.
Maaps expects a pre-aligned fasta file. 
Maaps currently does not support gene-coordinate annotations.

### Output
```
ID	  Codon	NTPos   	AA	 AAPos	 RowType
Den4	TAA	   -	   !	  3388	 Stop_Codon
Den4	TAG	   -	   !	  3414	 Stop_Codon
Den4	TAA	   -	   !	  3422	 Stop_Codon
Den4	TAG	   -	   !	  3473	 Stop_Codon

2055	TGA    -	  !    121	Stop_Codon
2055	RAC	 1927 	N/D  643	Non_Synonymous
```

Line breaks indicate the start of annotations for a new sequence; sequences are reported in the order they appeared in the input file.
If a sequence ID does not appear in the output file, that means no annotations were found for that sequence.

### Features

Maaps reports the following annotations: 

  * Non-synonymous degenerate codons
  * Stop codons in the middle of the sequence
  * Gaps
  * Ns

### Planned Features
The following features are planned:

  * Gene coordinate support accepting any of the following:
      * Genbank ID, genbank file, TSV file
  * Report frame shifts
  * Integrate stop-codon-checking with gene coordinate information
  * Report if the sequence/gene does not end with a stop codon
  * Automatic alignment

The following development features are planned:
  * Regression tests
  * Travis Integration
  * hpc reporting
  * hlint reporitng
  * liquid tags
