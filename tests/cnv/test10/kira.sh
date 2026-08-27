set -euo pipefail
rm -rf out
kira-bt cnv -- -s Q -c C -f af.tab.gz -o out in.noaf.vcf.gz
for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf
