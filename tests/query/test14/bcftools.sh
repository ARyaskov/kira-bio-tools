bcftools query -H -f '[%CHROM %POS  %SAMPLE %DP %GT\n]' in.vcf > out.bcf.vcf
