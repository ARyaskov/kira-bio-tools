bcftools query -f '[%POS\t%SAMPLE\t%GQ\n]' -i 'N_PASS(GQ<20)==1' in.vcf > out.bcf.vcf
