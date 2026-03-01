bcftools reheader -h reheader.2.hdr -f reheader.fai in.vcf.gz | bcftools view --no-version > out.bcf.vcf
