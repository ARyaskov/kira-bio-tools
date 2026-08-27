kira-bt query -f '%CHROM:%POS\t%N_PASS(GT="alt" & GQ>110)\t[\t%GT]\t[\t%GQ]\n' in.vcf > out.kira.vcf
