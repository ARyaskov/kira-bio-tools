bcftools convert --no-version --vcf-ids --hapsample2vcf in.hap,in.sample | grep -v '^##' > out.bcf.vcf
