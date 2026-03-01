bcftools gtcheck -e 0 in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.bcf.vcf
