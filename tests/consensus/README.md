# Consensus Test Suite

This suite validates kira-bt consensus behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt consensus
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s A > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s B > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s A -a N > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s B -a N > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - -m mask.bed > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -H 1 -m mask.bed > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -I -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - -I -m mask.bed --mask-with X -m mask.bed --mask-with lc > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -H 1 > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - -a . > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - -a . -i 'type="snp" || type="ref"' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - -a . -e 'MinDP>15' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - -a . -i 'ALT!="<DEL>"' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - -a . -e 'MinDP<15' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s smpl > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test17
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s smpl -a N > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test18
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test19
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - --mark-del - --mark-ins uc --mark-snv uc > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test20
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test21
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test22
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -H I --mark-ins lc --mark-snv lc > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test23
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test24
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test25
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -M . -S samples.txt > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test26
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -M . -s b > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test27
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -M . -s a > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test28
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test29
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - -M N -m mask.bed > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test30
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test31
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test32
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s - -m mask.bed -c out.chain > /dev/null; cat out.chain > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test33
- Command: kira-bt consensus -- in.bcf -f ref.fa -s - > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test34
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s S -H A > out.kira.vcf
- Checks: `-H A` mode chooses ALT in heterozygous `0/1`, keeps ALT in homozygous ALT, and does not alter homozygous REF sites.

test35
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s S -H LA > out.kira.vcf
- Checks: `-H LA` mode chooses the longer allele; when REF and ALT have equal length in het, tie is resolved to ALT as in bcftools.

test36
- Command: kira-bt consensus -- in.vcf.gz -f ref.fa -s S -H 1pIu > out.kira.vcf
- Checks: `NpIu` hybrid mode uses the Nth allele for phased genotypes and IUPAC ambiguity code for unphased genotypes.

## Pass Criteria

A test passes if:
1. kira.sh runs without errors.
2. out.kira.vcf matches out.kira.ref.vcf.

## Updating References

1. Rebuild kira-bt.
2. Run kira.sh in the target testN directory.
3. If behavior changes are expected, update out.kira.ref.vcf.
4. If bcftools.sh exists, update out.bcf.ref.vcf as upstream control.
