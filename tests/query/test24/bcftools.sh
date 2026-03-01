bcftools query -f '%POS\t%II[\t%FI]\n' -i 'sum(II)==6' in.vcf > out.bcf.vcf
