bcftools gtcheck in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v 'Time' > out.bcf.vcf
