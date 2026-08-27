bcftools gtcheck -t 11:33 -p A,D,A,E,D,E -u GT -e 10 in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.bcf.vcf
