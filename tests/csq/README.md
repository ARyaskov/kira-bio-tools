# Csq Test Suite

This suite validates kira-bt csq behavior and regression stability.

## Test Layout

In each testN directory:
- scenario input files
- kira.sh command for kira-bt csq
- optional bcftools.sh command for upstream reference generation
- out.kira.ref.vcf expected output used by CI
- out.kira.vcf latest produced output

## Test Case Matrix

test1
- Command: kira-bt csq -- -f ref.fa -g ann.gff in.vcf | gawk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {m=split($8,a,";"); for(i=1;i<=m;i++){ if(a[i] ~ /^(EXP|BCSQ)=/){ split(a[i],kv,"="); n=split(kv[2],v,","); asort(v); s=v[1]; for(j=2;j<=n;j++) s=s","v[j]; a[i]=kv[1]"="s } } info=a[1]; for(i=2;i<=m;i++) info=info";"a[i]; $8=info; print }' | bcftools query -f '%POS\t%REF\t%ALT\t%EXP\n%POS\t%REF\t%ALT\t%BCSQ\n\n' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test2
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | gawk 'BEGIN{FS=OFS="\t"} /^#/ {print; next} {m=split($8,a,";"); for(i=1;i<=m;i++){ if(a[i] ~ /^(EXP|BCSQ)=/){ split(a[i],kv,"="); n=split(kv[2],v,","); asort(v); s=v[1]; for(j=2;j<=n;j++) s=s","v[j]; a[i]=kv[1]"="s } } info=a[1]; for(i=2;i<=m;i++) info=info";"a[i]; $8=info; print }' | bcftools query -f '%POS\t%REF\t%ALT\t%EXP\n%POS\t%REF\t%ALT\t%BCSQ\n\n' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test3
- Command: kira-bt csq -- -f ref.fa -g ann.gff in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test4
- Command: kira-bt csq -- -f ref.fa -g ann.gff --ncsq 64 in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test5
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test6
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test7
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test8
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test9
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test10
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test11
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test12
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test13
- Command: kira-bt csq -- -f ref.fa -g ann.gff --ncsq 32 -p a in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test14
- Command: kira-bt csq -- -f ref.fa -g ann.gff in.vcf | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT:%BCSQ]\n' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test15
- Command: kira-bt csq -- -f ref.fa -g ann.gff --ncsq 16 -p r in.vcf | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT:%BCSQ]\n' > out.kira.vcf 
- Checks: scenario-specific behavior and stable output for this command.

test16
- Command: kira-bt csq -- -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT:%BCSQ]\n' > out.kira.vcf 
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
