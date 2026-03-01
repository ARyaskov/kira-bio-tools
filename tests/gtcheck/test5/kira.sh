kira-bt gtcheck -- -e 0 -u GT,GT in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf
