bgzip -c annots4.tab > annots4.tab.gz
tabix -f -s1 -b2 -e2 annots4.tab.gz
bcftools annotate -a annots4.tab.gz -c CHROM,POS,REF,ALT,+FA,+FR,+IA,+IR,+SA,+SR -o annotate4.bcf.vcf annotate4.vcf
