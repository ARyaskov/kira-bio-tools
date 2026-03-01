set -e
bcftools index --tbi -f in.vcf.gz; bcftools index -s in.vcf.gz > out.bcf.vcf
