kira-bt gtcheck -- -e 0 -P pairs.txt in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf
