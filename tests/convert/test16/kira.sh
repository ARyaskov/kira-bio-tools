kira-bt convert -- --no-version --gvcf2vcf -i 'FILTER="PASS"' -f ref.fa in.vcf.gz | grep -v '^##bcftools' > out.kira.vcf
