set -e
if bcftools polysomy 2>&1 | grep -q "Usage:   bcftools polysomy"; then
  bcftools polysomy -i in.vcf.gz > out.bcf.vcf
elif bcftools +polysomy -h >/dev/null 2>&1; then
  bcftools +polysomy in.vcf.gz -- -i > out.bcf.vcf
else
  echo "SKIP_UNSUPPORTED_POLYSOMY" > out.bcf.vcf
fi
