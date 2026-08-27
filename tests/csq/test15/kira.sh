kira-bt csq -- -f ref.fa -g ann.gff --ncsq 16 -p r in.vcf | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT:%BCSQ]\n' > out.kira.vcf
