set -e
bcftools index -f in.bcf; bcftools index -n in.bcf.csi > out.bcf.vcf
