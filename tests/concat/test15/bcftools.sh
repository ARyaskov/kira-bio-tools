bcftools concat --no-version -G -a -D concat.5.a.vcf.gz concat.5.b.vcf.gz concat.5.c.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.bcf.vcf
