bcftools csq -f ref.fa -g ann.gff --ncsq 32 -p a in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.bcf.vcf
