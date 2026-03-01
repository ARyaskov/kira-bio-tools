set -euo pipefail
rm -rf out
bcftools cnv -s Q -c C --regions-overlap 2 -R regions.bed -o out in.vcf.gz
for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.bcf.vcf
