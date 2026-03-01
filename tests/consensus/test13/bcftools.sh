bcftools consensus in.vcf.gz -f ref.fa -s - -a . -e 'MinDP>15' > out.bcf.vcf
