bcftools query -f '%POS  %NUM_TAG\n' -i 'COUNT(INFO/NUM_TAG)=2' in.vcf > out.bcf.vcf
