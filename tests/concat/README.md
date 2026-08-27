# Concat Test Suite

This suite validates kira-bt concat behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt concat
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt concat -- --no-version concat.1.a.vcf.gz concat.1.b.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt concat -- --no-version concat.1.a.bcf concat.1.b.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt concat -- --no-version -G concat.1.a.vcf.gz concat.1.b.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt concat -- --no-version -a concat.2.a.vcf.gz concat.2.b.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt concat -- --no-version -a concat.2.a.bcf concat.2.b.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt concat -- --no-version -aD concat.2.a.vcf.gz concat.2.b.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt concat -- --no-version -aD concat.2.a.bcf concat.2.b.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt concat -- --no-version -l --ligate-warn concat.3.a.vcf.gz concat.3.b.vcf.gz concat.3.0.vcf.gz concat.3.c.vcf.gz concat.3.d.vcf.gz concat.3.e.vcf.gz concat.3.f.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt concat -- --no-version -l --ligate-warn concat.3.a.bcf concat.3.b.bcf concat.3.0.bcf concat.3.c.bcf concat.3.d.bcf concat.3.e.bcf concat.3.f.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt concat -- --no-version -l concat.4.a.vcf.gz concat.4.b.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt concat -- --no-version -l concat.4.a.bcf concat.4.b.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt concat -- --no-version -l --ligate-warn concat.5.a.vcf.gz concat.5.b.vcf.gz concat.5.c.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt concat -- --no-version -l --ligate-warn concat.5.a.bcf concat.5.b.bcf concat.5.c.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt concat -- --no-version -l --ligate-force concat.5.a.bcf concat.5.b.bcf concat.5.c.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt concat -- --no-version -G -a -D concat.5.a.vcf.gz concat.5.b.vcf.gz concat.5.c.vcf.gz | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt concat -- --no-version -G -a -D concat.5.a.bcf concat.5.b.bcf concat.5.c.bcf | bcftools view | grep -v '^##bcftools_' > out.kira.vcf 
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
