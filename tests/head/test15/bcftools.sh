bcftools view --no-version -Ob in.vcf.gz | bcftools head -s1 -h2 > out.bcf.vcf
