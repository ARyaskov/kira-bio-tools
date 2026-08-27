bcftools filter --no-version -i 'sum(AD[*]) > FORMAT/DP' in.vcf -o out.bcf.vcf
