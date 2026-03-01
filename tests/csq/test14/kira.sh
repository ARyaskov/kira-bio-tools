kira-bt csq -- -f ref.fa -g ann.gff in.vcf | bcftools query -f '%POS\t%REF\t%ALT\t%INFO/BCSQ[\t%GT:%BCSQ]\n' > out.kira.vcf
