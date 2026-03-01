bcftools concat --no-version  concat.1.a.vcf.gz concat.1.b.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.bcf.vcf
