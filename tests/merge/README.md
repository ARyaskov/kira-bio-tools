# Merge Test Suite

This suite validates kira-bt merge behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt merge
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt merge merge.a.vcf.gz merge.b.vcf.gz -- --force-samples > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt merge merge.a.vcf.gz merge.b.vcf.gz merge.c.vcf.gz -- --force-samples > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt merge merge.a.vcf.gz merge.b.vcf.gz merge.c.vcf.gz -- --force-samples -Fx > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt merge merge.a.vcf.gz merge.b.vcf.gz merge.c.vcf.gz -- --force-samples -0 > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt merge merge.2.a.vcf.gz merge.2.b.vcf.gz -- --force-samples -m none > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt merge merge.2.a.vcf.gz merge.2.b.vcf.gz -- --force-samples -m both > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt merge merge.2.a.vcf.gz merge.2.b.vcf.gz -- --force-samples -m all > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt merge merge.3.a.vcf.gz merge.3.b.vcf.gz -- --force-samples -i TR:sum,TA:sum,TG:sum > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt merge merge.4.a.vcf.gz merge.4.b.vcf.gz -- --force-samples -m id > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt merge merge.5.a.vcf.gz merge.5.b.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt merge merge.6.a.vcf.gz merge.6.b.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt merge merge.8.a.vcf.gz merge.8.b.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt merge merge.8.a.vcf.gz merge.8.b.vcf.gz -- -i AN:sum,AC:sum > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt merge merge.9.a.vcf.gz merge.9.b.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt merge merge.9.a.vcf.gz merge.9.b.vcf.gz -- -i AN:sum,AC:sum > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt merge merge.10.a.vcf.gz merge.10.b.vcf.gz -- -m none > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test17
- Command: kira-bt merge merge.10.a.vcf.gz merge.10.b.vcf.gz -- -m both > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test18
- Command: kira-bt merge merge.10.a.vcf.gz merge.10.b.vcf.gz -- -m snp-ins-del > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test19
- Command: kira-bt merge merge.join.a.vcf.gz merge.join.b.vcf.gz -- -i AF:join > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test20
- Command: kira-bt merge merge.symbolic.1.a.vcf.gz merge.symbolic.1.b.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test21
- Command: kira-bt merge merge.multiallelics.1.a.vcf.gz merge.multiallelics.1.b.vcf.gz -- --merge none > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test22
- Command: kira-bt merge merge.multiallelics.1.a.vcf.gz merge.multiallelics.1.b.vcf.gz -- --merge both > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test23
- Command: kira-bt merge merge.phased.1.a.vcf.gz merge.phased.1.b.vcf.gz > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test24
- Command: kira-bt merge merge.7.a.vcf.gz merge.7.b.vcf.gz -- --force-samples > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test25
- Command: kira-bt merge merge.12.a.vcf.gz merge.12.b.vcf.gz -- --merge none > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test26
- Command: kira-bt merge merge.12.a.vcf.gz merge.12.b.vcf.gz -- --merge both > out.kira.vcf 
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
