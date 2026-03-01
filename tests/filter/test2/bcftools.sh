bcftools filter --no-version -e 'QUAL==59.2 || (INDEL=0 & (FMT/GQ=25 | FMT/DP=10))' -s Modified -S . in.vcf -o out.bcf.vcf
