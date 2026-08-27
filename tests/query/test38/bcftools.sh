bcftools query -f '%CHROM\t%POS\t%INFO/STR\n' -i 'INFO/STR=@query.string.2.1.txt' in.vcf > out.bcf.vcf
