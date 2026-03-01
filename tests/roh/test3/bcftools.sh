bcftools roh in.vcf.gz -Or -G30 --AF-file roh.1.tab.gz --ignore-homref | grep -v '^#' > out.bcf.vcf
