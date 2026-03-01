kira-bt query -f '%CHROM:%POS[\t%SAMPLE=%GT]\n' -e 'GT="mis"' -s 1,3,0 in.vcf > out.kira.vcf
