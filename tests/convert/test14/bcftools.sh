bcftools convert --no-version --hapsample2vcf in.hap,in.sample | grep -v '^##' > out.bcf.vcf
