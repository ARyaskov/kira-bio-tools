bcftools consensus in.vcf.gz -f ref.fa -s - -a . -i 'ALT!="<DEL>"' > out.bcf.vcf
