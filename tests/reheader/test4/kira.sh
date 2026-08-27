kira-bt reheader in.vcf.gz -- -o out.tmp.bcf -h reheader.hdr
bcftools view --no-version out.tmp.bcf > out.kira.vcf
