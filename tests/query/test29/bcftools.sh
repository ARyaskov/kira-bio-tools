bcftools query -f '[%POS\t%SAMPLE\t%GT\n]' -i 'N_PASS(GT="alt")==1' in.vcf > out.bcf.vcf
