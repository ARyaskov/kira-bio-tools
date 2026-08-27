bcftools query -f '[%GT]\n' -i 'GT~"0/[1-9]" || GT~"[1-9]/0"' in.vcf > out.bcf.vcf
