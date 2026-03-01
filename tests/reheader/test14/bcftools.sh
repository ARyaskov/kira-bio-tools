cat in.vcf.gz | bcftools reheader -s reheader.samples2 | bcftools view --no-version > out.bcf.vcf
