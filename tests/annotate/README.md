# Annotate Test Suite

This suite validates kira-bt annotate behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt annotate
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt annotate -a db.tab -c CHROM,FROM,TO,T_STR -o out.kira.vcf in.vcf
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt annotate -a db.vcf -c ID,QUAL,FILTER,INFO,FMT in.vcf -o out.kira.vcf
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: missing kira.sh
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt annotate -a annots4.tab -c CHROM,POS,REF,ALT,+FA,+FR,+IA,+IR,+SA,+SR annotate4.vcf -o annotate4.kira.vcf
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt annotate -a db.vcf.gz -c +INFO in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt annotate -a db.vcf.gz -c ALT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt annotate -a db.vcf.gz -c +ALT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt annotate -a db.vcf.gz -c ID,QUAL,FILTER,INFO,FMT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt annotate -a db.vcf.gz -c AAA:=IINT,FMT/BBB:=FMT/FINT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt annotate -a db.vcf.gz -c INFO/FILTER:=FILTER in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt annotate -a db.vcf.gz -c FMT/GT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt annotate -a db.vcf.gz -c STR,ID,QUAL,FILTER in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt annotate -a db.vcf.gz -c FMT/newGT:=GT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt annotate -a db.vcf.gz -c FMT/GT:=newGT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt annotate -a db.vcf.gz -c FILTER,INFO/FILTER:=./FILTER in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt annotate -a db.vcf.gz -c +FMT/GT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test17
- Command: kira-bt annotate -a db.vcf.gz -c XX in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test18
- Command: kira-bt annotate -a db.vcf.gz -c FMT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test19
- Command: kira-bt annotate -a db.vcf.gz -c ID,QUAL,+FILTER,+INFO,FMT/GT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test20
- Command: kira-bt annotate -a db.vcf.gz -c ID,QUAL,+FILTER,+INFO in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test21
- Command: kira-bt annotate -a db.vcf.gz -c ID,QUAL,INFO,FMT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test22
- Command: kira-bt annotate -a db.vcf.gz -c ID,QUAL,FILTER,INFO in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test23
- Command: kira-bt annotate -a db.vcf.gz -c ID,QUAL,FILTER in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test24
- Command: kira-bt annotate -a db.vcf.gz -c FILTER,INFO/FILTER:=FILTER,INFO/INFO_FILTER:=INFO/FILTER in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test25
- Command: kira-bt annotate -a db.vcf.gz -c INFO/FILTER:=FILTER,INFO/INFO_FILTER:=INFO/FILTER in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test26
- Command: kira-bt annotate -a db.vcf.gz -c INFO/FILTER:=./FILTER,FILTER in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test27
- Command: kira-bt annotate -a db.vcf.gz -c FILTER,INFO/FILTER:=./FILTER in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test28
- Command: kira-bt annotate -a db.vcf.gz -c INFO/ID:=ID,INFO/INFO_ID:=INFO/ID,ID,=ID:=INFO/ID in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test29
- Command: kira-bt annotate -a db.vcf.gz -c FMT/newGT:=GT,FMT/GT:=GT in.vcf.gz -o out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test30
- Command: kira-bt annotate -a db.vcf.gz -c FMT/newGT:=GT,ID in.vcf.gz -o out.kira.vcf 
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
