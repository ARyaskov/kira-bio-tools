kira-bt reheader in.vcf.gz -- -o out.tmp.bcf -s reheader.samples2
bcftools view --no-version out.tmp.bcf > out.kira.vcf
