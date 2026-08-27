set -e
bcftools index -f in.bcf; bcftools index -s in.bcf > out.bcf.vcf
