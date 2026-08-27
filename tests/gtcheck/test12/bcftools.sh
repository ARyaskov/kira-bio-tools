bcftools gtcheck -e 0 -s qry:E,D,C in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.bcf.vcf
