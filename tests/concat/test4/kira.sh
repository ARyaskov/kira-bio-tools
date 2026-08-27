kira-bt concat -- --no-version -a concat.2.a.vcf.gz concat.2.b.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf
