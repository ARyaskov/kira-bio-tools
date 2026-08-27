set -e
bcftools index --csi -f in.vcf.gz; bcftools index -s in.vcf.gz > out.bcf.vcf
