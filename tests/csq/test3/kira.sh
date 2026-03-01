kira-bt csq -- -f ref.fa -g ann.gff  in.vcf | bcftools query -f '[%TBCSQ\n]' | tr '\t,' '\n\n' | sed '/^$/d' | sort > out.kira.vcf
