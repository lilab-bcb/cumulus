### Drop-Seq Bundle
#### Create resource files for the Drop-Seq pipeline
Available at https://portal.firecloud.org/#methods/regev/dropseq_bundle/

Drop-Seq bundle requires the following files:

1. Reference genome sequence (FASTA file) 
1. Gene annotations (GTF file). Genes must have annotations with feature type 'exon' (column 3) in the GTF file in order to be used in alignment. 
The attribute keys transcript_id and gene_id in column 9 of the GTF file are required.


#### Single Species 
Input one set of matched FASTA and GTF files typically obtained from Ensembl, NCBI, or UCSC.

#### Multiple Species
This is similar to the single species case, but note that the species order of the FASTA and GTF files must match.
For example:

fasta_file: ["gs://fc-xxx/hg19.fa", "gs://fc-xxx/mm10.fa"]

gtf_file: ["gs://fc-xxx/hg19.gtf", "gs://fc-xxx/mm10.gtf"]




