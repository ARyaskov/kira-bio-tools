bcftools csq -f ref.fa -g ann.gff --unify-chr-names chr,chr,chr in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.bcf.vcf
