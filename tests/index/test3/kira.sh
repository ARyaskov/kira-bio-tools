set -e
kira-bt index in.vcf.gz -- --tbi -f; kira-bt index in.vcf.gz -- -n > out.kira.vcf
