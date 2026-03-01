bcftools gtcheck -e 0 --n-matches 4 in.vcf.gz | grep -v '^#' | grep -v '^INFO' | sort > out.bcf.vcf
