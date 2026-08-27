bcftools gtcheck -e 0 -p B,C,B,D in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.bcf.vcf
