bcftools gtcheck -e 0 in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.bcf.vcf
