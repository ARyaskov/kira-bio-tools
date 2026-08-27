bcftools query -f '[%POS  %SAMPLE  %AD\n]' -i 'FMT/AD[:0] < FMT/AD[:1]' in.vcf > out.bcf.vcf
