bcftools mpileup --no-version -f ref.fa -a DP,DV,-AD -r17:100-600 --gvcf 0,2,5 mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.bcf.vcf
