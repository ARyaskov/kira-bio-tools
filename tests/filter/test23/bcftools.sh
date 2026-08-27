bcftools filter --no-version -S . -e 'MAX(FORMAT/AO[0:])==4' in.vcf -o out.bcf.vcf
