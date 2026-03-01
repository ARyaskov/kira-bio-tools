# Stats Test Suite

This suite validates kira-bt stats behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt stats
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt stats stats.a.vcf.gz stats.b.vcf.gz -- -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt stats stats.a.vcf.gz stats.b.vcf.gz -- -s B > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt stats stats.counts.vcf.gz -- -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt stats stats.counts.vcf.gz -- -s - -i 'type="snp"' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt stats stats.vaf.vcf.gz -- -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt stats stats.counts.vcf.gz -- -s - -e 'type="snp"' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt stats stats.counts.vcf.gz -- -s - -1 > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt stats stats.counts.vcf.gz -- -s - -r 1 > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt stats stats.counts.vcf.gz -- -s - -r 2 > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt stats stats.counts.vcf.gz -- -s - -f PASS,. > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt stats stats.a.vcf.gz stats.b.vcf.gz -- -s - -r 1 > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt stats stats.a.vcf.gz stats.b.vcf.gz -- -s - -r 2 > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt stats stats.vaf.vcf.gz -- -s - -i 'QUAL>10' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt stats stats.vaf.vcf.gz -- -s - -e 'QUAL>10' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt stats stats.vaf.vcf.gz -- -s - -v > out.kira.vcf 
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
