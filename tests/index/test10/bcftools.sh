set -e
bcftools index --csi -f -o custom.csi in.vcf.gz; [ -s custom.csi ] && echo OK > out.bcf.vcf
