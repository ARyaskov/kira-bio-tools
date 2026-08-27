set -e
kira-bt index in.vcf.gz -- --csi -f; kira-bt index in.vcf.gz.csi -- -n > out.kira.vcf
