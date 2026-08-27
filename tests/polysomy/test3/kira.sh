set -e
if bcftools polysomy 2>&1 | grep -q "Usage:   bcftools polysomy" || bcftools +polysomy -h >/dev/null 2>&1; then
  kira-bt polysomy in.vcf.gz -- -r 1:100174876-100318245 > out.kira.vcf
else
  echo "SKIP_UNSUPPORTED_POLYSOMY" > out.kira.vcf
fi
