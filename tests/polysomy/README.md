# Polysomy Test Suite

This suite validates kira-bt polysomy behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt polysomy
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: set -e; if bcftools polysomy 2>&1 | grep -q "Usage: bcftools polysomy" || bcftools +polysomy -h >/dev/null 2>&1; then; kira-bt polysomy in.vcf.gz -- > out.kira.vcf; else; echo "SKIP_UNSUPPORTED_POLYSOMY" > out.kira.vcf; fi 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: set -e; if bcftools polysomy 2>&1 | grep -q "Usage: bcftools polysomy" || bcftools +polysomy -h >/dev/null 2>&1; then; kira-bt polysomy in.vcf.gz -- -s sample > out.kira.vcf; else; echo "SKIP_UNSUPPORTED_POLYSOMY" > out.kira.vcf; fi 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: set -e; if bcftools polysomy 2>&1 | grep -q "Usage: bcftools polysomy" || bcftools +polysomy -h >/dev/null 2>&1; then; kira-bt polysomy in.vcf.gz -- -r 1:100174876-100318245 > out.kira.vcf; else; echo "SKIP_UNSUPPORTED_POLYSOMY" > out.kira.vcf; fi 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: set -e; if bcftools polysomy 2>&1 | grep -q "Usage: bcftools polysomy" || bcftools +polysomy -h >/dev/null 2>&1; then; kira-bt polysomy in.vcf.gz -- -m 0.2 > out.kira.vcf; else; echo "SKIP_UNSUPPORTED_POLYSOMY" > out.kira.vcf; fi 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: set -e; if bcftools polysomy 2>&1 | grep -q "Usage: bcftools polysomy" || bcftools +polysomy -h >/dev/null 2>&1; then; kira-bt polysomy in.vcf.gz -- -f 0.10 > out.kira.vcf; else; echo "SKIP_UNSUPPORTED_POLYSOMY" > out.kira.vcf; fi 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: set -e; if bcftools polysomy 2>&1 | grep -q "Usage: bcftools polysomy" || bcftools +polysomy -h >/dev/null 2>&1; then; kira-bt polysomy in.vcf.gz -- -b 0.35 -c 0.7 -p 0.6 > out.kira.vcf; else; echo "SKIP_UNSUPPORTED_POLYSOMY" > out.kira.vcf; fi 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: set -e; if bcftools polysomy 2>&1 | grep -q "Usage: bcftools polysomy" || bcftools +polysomy -h >/dev/null 2>&1; then; kira-bt polysomy in.vcf.gz -- -i > out.kira.vcf; else; echo "SKIP_UNSUPPORTED_POLYSOMY" > out.kira.vcf; fi 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: set -e; if bcftools polysomy 2>&1 | grep -q "Usage: bcftools polysomy" || bcftools +polysomy -h >/dev/null 2>&1; then; kira-bt polysomy in.vcf.gz -- -v 2 > out.kira.vcf; else; echo "SKIP_UNSUPPORTED_POLYSOMY" > out.kira.vcf; fi 
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
