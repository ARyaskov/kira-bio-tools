bcftools concat --no-version  concat.1.a.bcf concat.1.b.bcf | bcftools view | grep -v '^##bcftools_' > out.bcf.vcf
