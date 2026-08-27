bcftools mpileup --no-version -f ref.fa -a -AD -t 1:100 --skip-all-set PAIRED,PROPER_PAIR,MREVERSE mpileup-filter.sam -o out.bcf.vcf
