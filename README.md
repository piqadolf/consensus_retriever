# consensus_retriever

## consensus_retriever is a small CLI application designed to get consensus DNA sequences from genome combined with Variant Call Format (VCF) file.
Pure python3 with only built-in packages.

### Usage
It takes genomic data in .fasta format and a corresponding .vcf file as an input.
The program has two modes. The first one incorporates variants in random places of random chromosomes.
The second mode allows user to specify which chromosomes and positions to investigate.

### Examples:
1) Random mode.
python3 consensus_retriever.py -g Felis_catus.Felis_catus_9.0.dna.toplevel.fa
	--vcf Contig.intervals_SNP.vcf -o random_cons.fa -f 50 -l 74 -n 200

2) Specified mode. The only difference from the previous mode is that you need to use --chr_pos flag.
You need to provide desired chromosomes and positions in the following format - 'chr1:123,456;chr2:789,890' (without quotes).

### Output:
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s1;allele=1
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACT
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s1;allele=2
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACT
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s2;allele=1
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACT
>chr=AANG04003642.1;pos=81454-81654;samples=fca.s2;allele=2
CCCTTTGATGTACCTGCAGTTCTGGCACAAATCCCAGCTGCAGACAGTCAGCTGGACTTCTAACCCTGCCCACT

Once you run the program, you have an sqlite database. Unless you delete it, all the subsequent runs will skip the genome parsing step and save you a substantial amount of time.
