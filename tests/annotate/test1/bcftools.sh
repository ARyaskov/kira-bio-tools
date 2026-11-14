bgzip -c db.tab > db.tab.gz
tabix -s1 -b2 -e3 db.tab.gz
printf '%s\n' '##INFO=<ID=T_STR,Number=.,Type=String,Description="Annotation field">' > annotation.hdr
bcftools annotate -a db.tab.gz -h annotation.hdr -c CHROM,FROM,TO,T_STR -o out.bcf.vcf in.vcf
