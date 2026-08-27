bcftools mpileup --no-version -f ref.fa -a DP,DPR,DV,DP4,INFO/DPR,SP,-AD -r17:100-600 mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.bcf.vcf
