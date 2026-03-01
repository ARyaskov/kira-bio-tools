set -e
bcftools index --tbi -f in.vcf.gz; bcftools index -n in.vcf.gz > out.bcf.vcf
