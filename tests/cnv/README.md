# Cnv Test Suite

This suite validates kira-bt cnv behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt cnv
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -r 1 -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -r 2:100-300 -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -R regions.bed -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -t 2:100-400 -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -T targets.bed -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C --regions-overlap 0 -R regions.bed -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C --regions-overlap 2 -R regions.bed -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -f af.tab.gz -o out in.noaf.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -a 0.7,0.8 -P 0.7 -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -b 0.9 -d 0.05,0.06 -e 1e-5 -k 0.3,0.25 -l 0.5 -L 5 -x 1e-8 -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: set -euo pipefail; rm -rf out; kira-bt cnv -- -s Q -c C -O 0.8 -o out in.vcf.gz; for f in summary.tab summary.Q.tab summary.C.tab cn.Q.tab cn.C.tab dat.Q.tab dat.C.tab; do if [ -f "out/$f" ]; then echo "###$f"; grep -v "^#" "out/$f"; fi; done > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

## Pass Criteria

A test passes if:
1. kira.sh runs without errors.
2. out.kira.vcf matches out.kira.ref.vcf.

## Updating References

1. Rebuild kira-bt.
2. Run kira.sh in the target testN directory.
3. If behavior changes are expected, update out.kira.ref.vcf.
4. If bcftools.sh exists, update out.bcf.ref.vcf as upstream control.
