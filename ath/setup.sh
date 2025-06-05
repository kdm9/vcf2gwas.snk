wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz
wget https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz.md5
md5sum -c 1001genomes_snp-short-indel_only_ACGTN.vcf.gz.md5

tabix 1001genomes_snp-short-indel_only_ACGTN.vcf.gz 
bcftools view --no-version \
    --samples-file samplelist.txt \
    -T regions_2m-5m-allchr.bed \
    -Ou \
    1001genomes_snp-short-indel_only_ACGTN.vcf.gz | \
bcftools annotate --no-version \
    -x 'ID,^INFO,^FORMAT/GT' -Ou | \
bcftools view  --no-version \
    -i 'MAC>5 & F_MISSING < 0.2' \
    --write-index \
    -Oz -o ath_filt-MAC5-MISS20.vcf.gz
