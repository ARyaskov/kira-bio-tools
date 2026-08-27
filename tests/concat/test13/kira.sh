kira-bt concat -- --no-version -l --ligate-warn concat.5.a.bcf concat.5.b.bcf concat.5.c.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf
