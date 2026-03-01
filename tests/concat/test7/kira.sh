kira-bt concat -- --no-version -aD concat.2.a.bcf concat.2.b.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf
