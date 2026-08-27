kira-bt concat -- --no-version -l concat.4.a.vcf.gz concat.4.b.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf
