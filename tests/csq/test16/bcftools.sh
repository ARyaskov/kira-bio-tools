bcftools csq -f ref.fa -g ann.gff --unify-chr-names 1 in.vcf | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT:%BCSQ]\n' > out.bcf.vcf
