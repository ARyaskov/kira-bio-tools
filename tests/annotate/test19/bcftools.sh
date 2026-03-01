bcftools annotate -a db.vcf.gz -c ID,QUAL,+FILTER,+INFO,FMT/GT  -o out.bcf.vcf in.vcf.gz
