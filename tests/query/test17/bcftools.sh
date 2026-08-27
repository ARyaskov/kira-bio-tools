bcftools query -f '%POS\t%REF\t%ALT\n' -i 'type!="snp"' in.vcf > out.bcf.vcf
