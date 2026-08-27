bcftools filter --no-version -i 'FORMAT/DP < sum(AD[*])' in.vcf -o out.bcf.vcf
