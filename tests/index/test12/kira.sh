set -e
kira-bt index in.bcf -- -f; kira-bt index in.bcf.csi -- -n > out.kira.vcf
