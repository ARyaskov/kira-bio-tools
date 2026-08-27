# Norm Test Suite

This suite validates kira-bt norm behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt norm
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt norm -f ref.fa in.vcf -o out.kira.vcf 
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
