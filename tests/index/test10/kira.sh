set -e
kira-bt index in.vcf.gz -- --csi -f -o custom.csi; [ -s custom.csi ] && echo OK > out.kira.vcf
