bcftools gtcheck -e 0 -P pairs.txt --distinctive-sites 3 in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.bcf.vcf
