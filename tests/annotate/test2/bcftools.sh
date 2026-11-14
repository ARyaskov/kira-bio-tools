bgzip -c db.vcf > db.vcf.gz
tabix -p vcf db.vcf.gz
bgzip -c in.vcf > in.vcf.gz
tabix -p vcf in.vcf.gz
bcftools annotate -a db.vcf.gz -c ID,QUAL,FILTER,INFO,FMT in.vcf.gz -o out.bcf.vcf
