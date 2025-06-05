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

## simulate traits

plink1.9 --make-bed --vcf ath_filt-MAC5-MISS20.vcf.gz --out ath_filt-MAC5-MISS20
alias simu="../sim/simu_linux"
simu --bfile ath_filt-MAC5-MISS20 --qt --causal-n 3 --hsq 0.5 --out ath_sim_h50
simu --bfile ath_filt-MAC5-MISS20 --qt --causal-n 3 --hsq 0.25 --out ath_sim_h25
simu --bfile ath_filt-MAC5-MISS20 --qt --causal-n 3 --hsq 0.1 --out ath_sim_h10
simu --bfile ath_filt-MAC5-MISS20 --qt --causal-n 3 --hsq 0.01 --out ath_sim_h01
 
csvtk join -f 1,2 -tT \
    <(csvtk rename -tT -f trait1 -n ath_sim_h50   ath_sim_h50.pheno) \
    <(csvtk rename -tT -f trait1 -n ath_sim_h25   ath_sim_h25.pheno) \
    <(csvtk rename -tT -f trait1 -n ath_sim_h10   ath_sim_h10.pheno) \
    <(csvtk rename -tT -f trait1 -n ath_sim_h01   ath_sim_h01.pheno) \
| csvtk cut  -tT  -f -2 \
| csvtk rename -tT -f FID -n indiv > simulated_phenos.tsv


