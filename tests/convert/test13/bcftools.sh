bcftools convert --no-version -G in.gen,in.sample | grep -v '^##' > out.bcf.vcf
