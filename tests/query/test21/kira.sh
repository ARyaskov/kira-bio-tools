kira-bt query -f '[%POS\t%SAMPLE\t%GT\t%AD\n]' -i 'GT="het" & binom(FMT/AD)>0.01' in.vcf > out.kira.vcf
