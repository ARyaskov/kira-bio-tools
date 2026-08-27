bcftools query -f '%POS\t%AD\n' -i 'binom(INFO/AD[0],INFO/AD[1])>0.01' in.vcf > out.bcf.vcf
