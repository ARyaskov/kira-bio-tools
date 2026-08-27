bcftools mpileup --no-version -f ref.fa -a DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -r17:100-600 mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.bcf.vcf
