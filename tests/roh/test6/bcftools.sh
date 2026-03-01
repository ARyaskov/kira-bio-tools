bcftools roh in.vcf.gz -G30 --AF-dflt 0.4 -r 1:100174876-100318245 --ignore-homref --include-noalt | grep -v '^#' > out.bcf.vcf
