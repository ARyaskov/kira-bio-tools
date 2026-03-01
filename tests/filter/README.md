# Filter Test Suite

This suite validates kira-bt filter behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt filter
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt filter --no-version --soft-filter XX -m x -g2 -G2 in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt filter --no-version -e 'QUAL==59.2 || (INDEL=0 & (FMT/GQ=25 | FMT/DP=10))' -s Modified -S . in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt filter --no-version -e 'INFO/DP=19' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt filter --no-version -e 'INFO/DP=19' -s XX in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt filter --no-version -e 'INFO/DP=19' -s XX -m + in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt filter --no-version -e 'INFO/DP=19' -s XX -m x in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt filter --no-version -e 'INFO/DP=19' -s XX -m +x in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt filter --no-version -e 'FMT/GT="0/2"' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt filter --no-version -i 'FMT/GT="0/0" && AC[*]=2' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt filter --no-version -i 'AC[*]=2 && FMT/GT="0/0"' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt filter --no-version -i 'ALT="."' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt filter --no-version -S . -i 'FORMAT/TEST3<25' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt filter --no-version -S . -i 'FORMAT/TEST4<25' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt filter --no-version -i 'GT="HOM"' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt filter --no-version -i 'GT="HET"' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt filter --no-version -i 'GT="HAP"' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test17
- Command: kira-bt filter --no-version -i 'AD[:1]=11' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test18
- Command: kira-bt filter --no-version -i 'AD[1:]=11' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test19
- Command: kira-bt filter --no-version -i 'FR[0:1]=11' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test20
- Command: kira-bt filter --no-version -i 'F_MISSING>=0.2' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test21
- Command: kira-bt filter --no-version -i 'F_PASS(GT=="mis")>=0.2' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test22
- Command: kira-bt filter --no-version -m x -s + -g2:mnp,indel,other in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test23
- Command: kira-bt filter --no-version -S . -e 'MAX(FORMAT/AO[0:])==4' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test24
- Command: kira-bt filter --no-version -S . -e 'SMPL_MAX(FORMAT/AO)==4' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test25
- Command: kira-bt filter --no-version -S . -e 'sMIN(FORMAT/AO)==2' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test26
- Command: kira-bt filter --no-version -S . -e 'sSUM(FORMAT/AO)==11' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test27
- Command: kira-bt filter --no-version -i 'QUAL/FMT/AD==55' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test28
- Command: kira-bt filter --no-version -i 'QUAL/INFO/AD==10' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test29
- Command: kira-bt filter --no-version -i 'sum(AD[*]) > FORMAT/DP' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test30
- Command: kira-bt filter --no-version -i 'FORMAT/DP < sum(AD[*])' in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test31
- Command: kira-bt filter --no-version --soft-filter xxx --mask 2:1005-1008 --mask-overlap 0 in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test32
- Command: kira-bt filter --no-version --soft-filter xxx --mask 2:1005-1008 --mask-overlap 1 in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test33
- Command: kira-bt filter --no-version --soft-filter xxx --mask 2:1005-1008 --mask-overlap 2 in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test34
- Command: kira-bt filter --no-version --soft-filter xxx --mask ^2:1005-1008 in.vcf -o out.kira.vcf 
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
