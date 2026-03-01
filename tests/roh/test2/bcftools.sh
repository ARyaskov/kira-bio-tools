bcftools roh in.vcf.gz -Or -G30 --AF-file roh.1.tab.gz | grep -v '^#' > out.bcf.vcf
