bcftools reheader -o out.tmp.bcf -h reheader.hdr in.vcf.gz
bcftools view --no-version out.tmp.bcf > out.bcf.vcf
