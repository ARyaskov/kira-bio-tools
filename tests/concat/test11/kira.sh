kira-bt concat -- --no-version -l concat.4.a.bcf concat.4.b.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf
