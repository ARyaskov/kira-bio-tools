bcftools query -f '[%SAMPLE %DP\n]' -i 'DP=1 || DP=2' in.vcf > out.bcf.vcf
