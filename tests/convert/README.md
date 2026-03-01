# Convert Test Suite

This suite validates kira-bt convert behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt convert
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt convert -- -g -,. in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt convert -- -g -,. --vcf-ids in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt convert -- -g -,. --vcf-ids --3N6 in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt convert -- -g .,- in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt convert -- -g -,. --tag PL in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt convert -- -h -,.,. in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt convert -- -h .,-,. --vcf-ids in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt convert -- -h .,.,- in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt convert -- --hapsample -,. in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt convert -- --hapsample -,. --vcf-ids in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt convert -- --hapsample .,- in.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt convert -- --no-version --vcf-ids -G in.gen,in.sample | grep -v '^##' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt convert -- --no-version -G in.gen,in.sample | grep -v '^##' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt convert -- --no-version --hapsample2vcf in.hap,in.sample | grep -v '^##' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt convert -- --no-version --vcf-ids --hapsample2vcf in.hap,in.sample | grep -v '^##' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt convert -- --no-version --gvcf2vcf -i 'FILTER="PASS"' -f ref.fa in.vcf.gz | grep -v '^##bcftools' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test17
- Command: kira-bt convert -- --no-version -c ID,CHROM,POS,AA -s SAMPLE1 -f ref.fa --tsv2vcf in.tsv > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test18
- Command: kira-bt convert -- --no-version -c -,CHROM,POS,REF,ALT -f ref.fa --tsv2vcf in.tsv > out.kira.vcf 
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
