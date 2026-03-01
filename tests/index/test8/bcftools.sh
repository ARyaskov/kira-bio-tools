set -e
bcftools index -f in.bcf; bcftools index -n in.bcf > out.bcf.vcf
