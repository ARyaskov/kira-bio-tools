bcftools query -f '%POS\t%REF\t%ALT\t%ILEN\n' -i 'ILEN==1' in.vcf > out.bcf.vcf
