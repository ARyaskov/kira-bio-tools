# Gtcheck Test Suite

This suite validates kira-bt gtcheck behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt gtcheck
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt gtcheck -- in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt gtcheck -- in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt gtcheck -- -E 0 in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt gtcheck -- -e 0 in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt gtcheck -- -e 0 -u GT,GT in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt gtcheck -- -e 0 -u PL,PL in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt gtcheck -- -e 0 -p s1,s1 in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt gtcheck -- -e 0 in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt gtcheck -- -e 0 in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt gtcheck -- -e 0 -P pairs.txt in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt gtcheck -- -e 0 --n-matches 4 in.vcf.gz | grep -v '^#' | grep -v '^INFO' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt gtcheck -- -e 0 -s qry:E,D,C in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt gtcheck -- -e 0 -s qry:B -s gt:D in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt gtcheck -- -e 0 -s qry:B -s gt:D,C in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt gtcheck -- -e 0 -p B,C,B,D in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt gtcheck -- -e 0 -u GT,GT -H in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test17
- Command: kira-bt gtcheck -- -e 0 -P pairs.txt --distinctive-sites 3 in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test18
- Command: kira-bt gtcheck -- in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test19
- Command: kira-bt gtcheck -- --n-matches 2 in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test20
- Command: kira-bt gtcheck -- in.vcf.gz -g gts.vcf.gz | grep -v '^#' | grep -v 'Time' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test21
- Command: kira-bt gtcheck -- -p A,B,B,C in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test22
- Command: kira-bt gtcheck -- -t 11:33 -p A,D,A,E,D,E -u GT -e 10 in.vcf.gz | grep -v '^#' | grep -v '^INFO' > out.kira.vcf 
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
