# Mpileup Test Suite

This suite validates kira-bt mpileup behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt mpileup
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt mpileup iupac.bam -o out.kira.vcf -- -f ref.fa -r 11:10-20 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -r17:100-102,17:102-103,17:103-104,17:104-105,17:100-105 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -r17:100-150 -a -AD 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a DP,DV -r17:100-600 -a -AD 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt mpileup mpileup.1.bam -o out.kira.vcf -- -f ref.fa -B --ff 0x14 -r17:1050-1060 -a -AD 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a DP,DPR,DV,DP4,INFO/DPR,SP,-AD -r17:100-600 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR -r17:100-600 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a DP,DV,-AD -r17:100-600 --gvcf 0,2,5 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a -AD -r17:100-150 -s HG00101,HG00102 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a -AD -r17:100-150 -S mplp.samples 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a -AD -r17:100-150 -s ^HG00101,HG00102 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a -AD -t17:100-150 -S mplp.9.samples 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a -AD -t17:100-150 -G mplp.10.samples 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt mpileup mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a -AD 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt mpileup mpileup.3.bam mpileup.4.bam -o out.kira.vcf -- -f ref.fa -a -AD -s HG00102 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt mpileup indel-AD.1.bam -o out.kira.vcf -- -f ref.fa -a AD 
- Checks: scenario-specific behavior and stable output for this command.

test17
- Command: kira-bt mpileup indel-AD.2.bam -o out.kira.vcf -- -f ref.fa -a AD -r 11:75 
- Checks: scenario-specific behavior and stable output for this command.

test18
- Command: kira-bt mpileup indel-AD.2.bam -o out.kira.vcf -- -f ref.fa -a AD -r 11:75 --ambig-reads incAD 
- Checks: scenario-specific behavior and stable output for this command.

test19
- Command: kira-bt mpileup indel-AD.2.bam -o out.kira.vcf -- -f ref.fa -a AD -r 11:75 --ambig-reads incAD0 
- Checks: scenario-specific behavior and stable output for this command.

test20
- Command: kira-bt mpileup mpileup-SCR.bam -o out.kira.vcf -- -f ref.fa -a -AD,INFO/SCR,FMT/SCR 
- Checks: scenario-specific behavior and stable output for this command.

test21
- Command: kira-bt mpileup mpileup-filter.sam -o out.kira.vcf -- -f ref.fa -a -AD -t 1:100 --skip-all-set PAIRED,PROPER_PAIR,MREVERSE 
- Checks: scenario-specific behavior and stable output for this command.

test22
- Command: kira-bt mpileup mpileup-filter.sam -o out.kira.vcf -- -f ref.fa -a -AD -t 1:100 --skip-any-unset READ1 
- Checks: scenario-specific behavior and stable output for this command.

test23
- Command: kira-bt mpileup annot-NMBZ.1.bam -o out.kira.vcf -- -f ref.fa -a -AD,INFO/NMBZ -r chr19:69-99 
- Checks: scenario-specific behavior and stable output for this command.

test24
- Command: kira-bt mpileup annot-NMBZ.2.bam -o out.kira.vcf -- -f ref.fa -a -AD,INFO/NMBZ -r chr6:75 
- Checks: scenario-specific behavior and stable output for this command.

test25
- Command: kira-bt mpileup annot-NMBZ.3.1.bam annot-NMBZ.3.2.bam -o out.kira.vcf -- -f ref.fa -a -AD,INFO/NMBZ -r chr16:75 
- Checks: scenario-specific behavior and stable output for this command.

test26
- Command: kira-bt mpileup mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a GQ,PL 
- Checks: scenario-specific behavior and stable output for this command.

test27
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -r17:100-150 -a GQ,PL,-AD 
- Checks: scenario-specific behavior and stable output for this command.

test28
- Command: kira-bt mpileup mpileup.1.bam mpileup.2.bam mpileup.3.bam -o out.kira.vcf -- -f ref.fa -a GQ,PL,-AD -r17:100-150 -s HG00101,HG00102 
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
