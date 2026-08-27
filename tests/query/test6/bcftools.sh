bcftools query -f '%POS[ %GL]\n' -i 'min(abs(GL[*:0]))=10' in.vcf > out.bcf.vcf
