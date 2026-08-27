bcftools query -f '%CHROM\t%POS\t%INFO\t%FORMAT\n' -s D,C in.vcf > out.bcf.vcf
