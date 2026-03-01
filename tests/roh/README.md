# Roh Test Suite

This suite validates kira-bt roh behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt roh
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt roh in.vcf.gz -- -Or -G30 --AF-dflt 0.4 | grep -v '^#' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt roh in.vcf.gz -- -Or -G30 --AF-file roh.1.tab.gz | grep -v '^#' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt roh in.vcf.gz -- -Or -G30 --AF-file roh.1.tab.gz --ignore-homref | grep -v '^#' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt roh in.vcf.gz -- -G30 --AF-dflt 0.4 -r 1:100174876-100318245 | grep -v '^#' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt roh in.vcf.gz -- -G30 --AF-dflt 0.4 -r 1:100174876-100318245 --ignore-homref | grep -v '^#' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt roh in.vcf.gz -- -G30 --AF-dflt 0.4 -r 1:100174876-100318245 --ignore-homref --include-noalt | grep -v '^#' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt roh in.vcf.gz -- -G30 --AF-dflt 0.4 -r 1:100174876-100318245 --include-noalt | grep -v '^#' > out.kira.vcf 
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
