cat in.vcf.gz | kira-bt reheader - -- -s reheader.samples2 | bcftools view --no-version > out.kira.vcf
