bcftools query -f '%POS[\t%GT]\n' -i 'COUNT(GT="het")=1' in.vcf > out.bcf.vcf
