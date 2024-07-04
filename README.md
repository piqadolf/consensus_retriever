# consensus_retriever

## consensus_retriever is a small CLI wrote in python3 with only built-in packages.

### Usage
It takes genomic data in .fasta format and a corresponding .vcf file as an input.

### Example:
python3 consensus_retriever.py -g Felis_catus.Felis_catus_9.0.dna.toplevel.fa
	--vcf Contig.intervals_SNP.vcf -o random_cons.fa -f 50 -l 75

### Output:
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s1;allele=1
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACTA
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s1;allele=2
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACTA
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s2;allele=1
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACTA
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s2;allele=2
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACTA

