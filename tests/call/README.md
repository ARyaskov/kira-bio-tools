# Call Test Suite

This suite validates kira-bt call behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt call
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt call in.vcf -o out.kira.vcf -- -cv 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt call in.vcf -o out.kira.vcf -- -m -s sample3 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt call in.vcf -o out.kira.vcf -- -m --ploidy-file call-ploidy.1.txt -S call-ploidy.1.ped 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mg0 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv -S mpileup.3.samples 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv -S mpileup.4.samples 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv -S mpileup.5.samples 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv --ploidy-file mpileup.ploidy -S mpileup.samples 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv --ploidy-file mpileup.ploidy -S mpileup.ped 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv --ploidy-file mpileup.ploidy -S mpileup.2.samples 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv -G mpileup.hwe.samples --group-samples-tag AD 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv -G - --group-samples-tag AD 
- Checks: scenario-specific behavior and stable output for this command.

test17
- Command: kira-bt call in.vcf -o out.kira.vcf -- -mv -F AN_POP,AC_POP 
- Checks: scenario-specific behavior and stable output for this command.

test18
- Command: kira-bt call in.vcf -o out.kira.vcf -- -m 
- Checks: scenario-specific behavior and stable output for this command.

test19
- Command: kira-bt call in.vcf -o out.kira.vcf -- -m -G call.af-fixation.txt 
- Checks: scenario-specific behavior and stable output for this command.

test20
- Command: kira-bt call in.vcf -o out.kira.vcf -- -m -G call.af-fixation.txt -a GP,GQ 
- Checks: scenario-specific behavior and stable output for this command.

test21
- Command: kira-bt call in.vcf -o out.kira.vcf -- -m -G call.af-fixation.txt -a GP,GQ -s NA07051_CEU,NA12843_CEU 
- Checks: scenario-specific behavior and stable output for this command.

test22
- Command: kira-bt call in.vcf -o out.kira.vcf -- -m -G call.af-fixation.txt -a GP,GQ -S samples.txt 
- Checks: scenario-specific behavior and stable output for this command.

test23
- Command: kira-bt call in.vcf -o out.kira.vcf -- -cv -a GP,GQ 
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
