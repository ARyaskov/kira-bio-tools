bcftools gtcheck in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.bcf.vcf
