bcftools consensus in.vcf.gz -f ref.fa -s - -m mask.bed -c out.chain > /dev/null; cat out.chain > out.bcf.vcf
