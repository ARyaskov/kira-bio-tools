set -e
bcftools index --csi -f in.vcf.gz; bcftools index -n in.vcf.gz.csi > out.bcf.vcf
