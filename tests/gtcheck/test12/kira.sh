kira-bt gtcheck -- -e 0 -s qry:E,D,C in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf
