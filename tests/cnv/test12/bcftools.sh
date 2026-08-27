set -euo pipefail
rm -rf out
bcftools cnv -s Q -c C -b 0.9 -d 0.05,0.06 -e 1e-5 -k 0.3,0.25 -l 0.5 -L 5 -x 1e-8 -o out in.vcf.gz
for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.bcf.vcf
