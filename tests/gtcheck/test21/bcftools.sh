bcftools gtcheck -p A,B,B,C in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.bcf.vcf
