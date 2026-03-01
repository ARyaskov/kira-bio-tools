bcftools annotate -a db.vcf.gz -c INFO/ID:=ID,INFO/INFO_ID:=INFO/ID,ID,=ID:=INFO/ID  -o out.bcf.vcf in.vcf.gz
