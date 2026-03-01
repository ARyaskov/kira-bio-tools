# Query Test Suite

This suite validates kira-bt query behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt query
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP4\t%AN[\t%GT\t%TGT]\n' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt query -f '%RSX\t%VKX\n' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt query -f '%POS %REF %ALT\n' -i 'REF~"C" && ALT[*]~"CT"' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt query -f '%POS %REF %ALT\n' -i 'N_ALT=2' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt query -f '%POS %AN\n' -i 'AN!=2*N_SAMPLES' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt query -f '%POS[ %GL]\n' -i 'min(abs(GL[*:0]))=10' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt query -f '%POS[ %GT]\n' -i 'AC[0]=3' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt query -f '%POS[ %GT]\n' -i 'MAC[0]=1' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt query -f '%LINE' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt query -f '%CHROM %POS\n' -i 'CHROM="4"' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt query -f '[%POS\t%SAMPLE\t%GQ\n]' -i 'N_PASS(GQ<20)==1' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt query -f '%CHROM\t%POS\t%INFO\t%FORMAT\n' -s D,C in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt query -f '%CHROM:%POS\t%N_PASS(GT="alt" & GQ>110)\t[\t%GT]\t[\t%GQ]\n' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt query -H -f '[%CHROM %POS %SAMPLE %DP %GT\n]' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt query -HH -f '%CHROM %POS[ %SAMPLE][ %DP][ %GT]' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt query -f '%POS\t%REF\t%ALT\n' -i 'type="snp"' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test17
- Command: kira-bt query -f '%POS\t%REF\t%ALT\n' -i 'type!="snp"' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test18
- Command: kira-bt query -f '%POS\t%REF\t%ALT\n' -i 'INFO/TYPE="xxx"' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test19
- Command: kira-bt query -f '[%GT]\n' -i 'GT~"0/[1-9]" || GT~"[1-9]/0"' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test20
- Command: kira-bt query -f '%POS[\t%GT]\n' -i 'COUNT(GT="het")=1' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test21
- Command: kira-bt query -f '[%POS\t%SAMPLE\t%GT\t%AD\n]' -i 'GT="het" & binom(FMT/AD)>0.01' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test22
- Command: kira-bt query -f '%POS\t%AD\n' -i 'binom(INFO/AD[0],INFO/AD[1])>0.01' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test23
- Command: kira-bt query -f '%POS\n' -i 'POS==16777217 || POS==33554432 || POS=118673904' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test24
- Command: kira-bt query -f '%POS\t%II[\t%FI]\n' -i 'sum(II)==6' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test25
- Command: kira-bt query -f '%POS\t%II[\t%FI]\n' -i 'median(FORMAT/FI)==1.5' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test26
- Command: kira-bt query -f '%POS\t%REF\t%ALT\t%ILEN\n' -i 'ILEN==1' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test27
- Command: kira-bt query -f '[%POS %SAMPLE %AD\n]' -i 'FMT/AD[:0] < FMT/AD[:1]' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test28
- Command: kira-bt query -f '%POS %NUM_TAG\n' -i 'COUNT(INFO/NUM_TAG)=2' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test29
- Command: kira-bt query -f '[%POS\t%SAMPLE\t%GT\n]' -i 'N_PASS(GT="alt")==1' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test30
- Command: kira-bt query -f '%CHROM:%POS[\t%SAMPLE=%GT]\n' -e 'GT="mis"' -s 1,3,0 in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test31
- Command: kira-bt query -f '%ID\n' -i 'ID~"s12"' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test32
- Command: kira-bt query -f '[%SAMPLE %GT\n]' -S query.smpl.11.txt in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test33
- Command: kira-bt query -l -S query.smpl.11.txt in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test34
- Command: kira-bt query -f '[%SAMPLE %GT\n]' -S ^query.smpl.11.txt in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test35
- Command: kira-bt query -l -S ^query.smpl.11.txt in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test36
- Command: kira-bt query -f '%CHROM\t%POS\t%CLNREVSTAT\n' -i 'CLNREVSTAT="criteria_provided,_conflicting_interpretations"' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test37
- Command: kira-bt query -f '%CHROM\t%POS\t%CLNREVSTAT\n' -i 'CLNREVSTAT="criteria_provided" && CLNREVSTAT="_conflicting_interpretations"' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test38
- Command: kira-bt query -f '%CHROM\t%POS\t%INFO/STR\n' -i 'INFO/STR=@query.string.2.1.txt' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test39
- Command: kira-bt query -f '[%SAMPLE %DP\n]' -i 'DP=1 || DP=2' in.vcf > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test40
- Command: kira-bt query -f '%CHROM %POS %ID %REF %ALT %QUAL %FILTER \t %INFO/CHROM %INFO/POS %INFO/ID %INFO/REF %INFO/ALT %INFO/QUAL %INFO/FILTER' in.vcf > out.kira.vcf 
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
