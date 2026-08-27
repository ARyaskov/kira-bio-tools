bcftools annotate -a db.vcf.gz -c INFO/FILTER:=FILTER,INFO/INFO_FILTER:=INFO/FILTER  -o out.bcf.vcf in.vcf.gz
