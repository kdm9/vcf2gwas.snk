export PATH=.opt/bin/:$PATH

rm -rf sim_1k1g* 1k1g_chr1_1m-2m* v2g.wd/
plink --make-bed --allow-extra-chr --biallelic-only --mind 1.0 --allow-no-sex --double-id --set-missing-var-ids '@:#' --vcf 1k1g_chr1:1000000-2000000_miss10.vcf.gz --out 1k1g_chr1_1m-2m

simu --bfile 1k1g_chr1_1m-2m --qt --causal-n 20 --hsq 0.6 --out sim_1k1g_chr1_h60
plink  --bfile 1k1g_chr1_1m-2m --out sim_1k1g_chr1_h60  --recode --make-bed --pheno sim_1k1g_chr1_h60.pheno
gemma -bfile sim_1k1g_chr1_h60 -gk 1 -o sim_1k1g_chr1_h60.gk -outdir .
gemma -bfile sim_1k1g_chr1_h60 -lmm 4 -miss 1.0 -k sim_1k1g_chr1_h60.gk.cXX.txt -o sim_1k1g_chr1_h60.gemma -outdir .
./plot_gemma.R --highlight-snps sim_1k1g_chr1_h60.1.causals sim_1k1g_chr1_h60.gemma.assoc.txt
csvtk rename -f IID,trait1 -n indiv,trait_h60 -tT sim_1k1g_chr1_h60.pheno | 
    csvtk cut -f indiv,trait_h60 -tT > sim_1k1g_chr1_h60.pheno.tsv

simu--bfile 1k1g_chr1_1m-2m --qt --causal-n 20 --hsq 0.95 --out sim_1k1g_chr1_h95
plink1.9  --bfile 1k1g_chr1_1m-2m --out sim_1k1g_chr1_h95  --recode --make-bed --pheno sim_1k1g_chr1_h95.pheno
gemma -bfile sim_1k1g_chr1_h95 -gk 1 -o sim_1k1g_chr1_h95.gk -outdir .
gemma -bfile sim_1k1g_chr1_h95 -lmm 4 -miss 1.0 -k sim_1k1g_chr1_h95.gk.cXX.txt -o sim_1k1g_chr1_h95.gemma -outdir .
./plot_gemma.R --highlight-snps sim_1k1g_chr1_h95.1.causals  sim_1k1g_chr1_h95.gemma.assoc.txt

csvtk rename -f IID,trait1 -n indiv,trait_h95 -tT sim_1k1g_chr1_h95.pheno | 
    csvtk cut -f indiv,trait_h95 -tT > sim_1k1g_chr1_h95.pheno.tsv
csvtk join -f indiv -tT sim_1k1g_chr1_h60.pheno.tsv sim_1k1g_chr1_h95.pheno.tsv  > simulated_phenotypes.tsv

snakemake --snakefile vcf2gwas.snk -j8

cat v2g.wd/out/trait_pves.tsv

./plot_gemma.R --highlight-snps sim_1k1g_chr1_h60.1.causals  v2g.wd/out/trait_h60.assoc.txt
./plot_gemma.R --highlight-snps sim_1k1g_chr1_h95.1.causals  v2g.wd/out/trait_h95.assoc.txt

gcta64 --bfile 1k1g_chr1_1m-2m --maf 0.005 --make-grm --out gcta_grm
gcta64 --bfile 1k1g_chr1_1m-2m --maf 0.005 --make-grm --out gcta_grm
gcta64 --grm gcta_grm --pheno sim_1k1g_chr1_h60.pheno --reml --out gcta_h60_greml --thread-num 12
gcta64 --grm gcta_grm --pheno sim_1k1g_chr1_h95.pheno --reml --out gcta_h95_greml --thread-num 12

gcta64 --grm gcta_grm --make-bK-sparse 0.001 --out gcta_grm_sparse001
gcta64 --grm-sparse gcta_grm_sparse001  --fastGWA-mlm  --pheno sim_1k1g_chr1_h60.pheno --out gcta_h60_fastGWA-mlm --thread-num 12
