# consensus_retriever

## consensus_retriever is a small CLI application designed to get consensus DNA sequences from genome combined with Variant Call Format (VCF) file.
Pure python3 with only built-in packages.

### Usage
It takes genomic data in .fasta format and a corresponding .vcf file as an input.
It takes random slices from genome sequences and substitutes SNP positions with correspponding alleles from .vcf file.

### Example:
python3 consensus_retriever.py -g Felis_catus.Felis_catus_9.0.dna.toplevel.fa
	--vcf Contig.intervals_SNP.vcf -o random_cons.fa -f 50 -l 74 -n 200

### Output:
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s1;allele=1
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACT
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s1;allele=2
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACT
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s2;allele=1
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACT
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s2;allele=2
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACT
