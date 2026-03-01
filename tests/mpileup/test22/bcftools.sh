bcftools mpileup --no-version -f ref.fa -a -AD -t 1:100 --skip-any-unset READ1 mpileup-filter.sam -o out.bcf.vcf
