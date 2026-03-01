set -e
if bcftools polysomy 2>&1 | grep -q "Usage:   bcftools polysomy"; then
  bcftools polysomy -b 0.35 -c 0.7 -p 0.6 in.vcf.gz > out.bcf.vcf
elif bcftools +polysomy -h >/dev/null 2>&1; then
  bcftools +polysomy in.vcf.gz -- -b 0.35 -c 0.7 -p 0.6 > out.bcf.vcf
else
  echo "SKIP_UNSUPPORTED_POLYSOMY" > out.bcf.vcf
fi
