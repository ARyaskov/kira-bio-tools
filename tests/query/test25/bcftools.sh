bcftools query -f '%POS\t%II[\t%FI]\n' -i 'median(FORMAT/FI)==1.5' in.vcf > out.bcf.vcf
