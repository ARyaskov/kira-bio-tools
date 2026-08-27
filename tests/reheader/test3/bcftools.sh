bcftools reheader -o out.tmp.bcf -s reheader.samples2 in.vcf.gz
bcftools view --no-version out.tmp.bcf > out.bcf.vcf
