bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP4\t%AN[\t%GT\t%TGT]\n' in.vcf > out.bcf.vcf
