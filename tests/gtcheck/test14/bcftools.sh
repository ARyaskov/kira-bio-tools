bcftools gtcheck -e 0 -s qry:B -s gt:D,C in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.bcf.vcf
