# Sort Test Suite

This suite validates kira-bt sort behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt sort
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt sort in.vcf -o out.kira.vcf -- -m 0 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt sort in.vcf -o out.kira.vcf -- -m 1000 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt sort in.vcf -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt sort in.vcf -o out.kira.vcf -- -m 1M 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt sort in.vcf -o out.kira.vcf -- -m 10K 
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
