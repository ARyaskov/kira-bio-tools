bcftools convert --no-version --vcf-ids -G in.gen,in.sample | grep -v '^##' > out.bcf.vcf
