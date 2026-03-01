bcftools query -f '%POS %REF %ALT\n' -i 'REF~"C" && ALT[*]~"CT"' in.vcf > out.bcf.vcf
