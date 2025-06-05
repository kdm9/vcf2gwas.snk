python3 fst_vcfsim.py  --fst ../test/eucalpts.pops.txt --pop-sizes 100,100,100,10 --num-variants 10000 --output-prefix euc_sim --verbose
plink1.9 --make-bed --vcf OUT_euc3.vcf --out euc_sim
./simu_linux --bfile euc_sim --qt --causal-n 10 --hsq 0.1 --out euc_sim_h10
./simu_linux --bfile euc_sim --qt --causal-n 10 --hsq 0.01 --out euc_sim_h01

plink1.9  --bfile euc_sim --out euc_sim_h10  --recode --make-bed --pheno euc_sim_h10.pheno
plink1.9  --bfile euc_sim --out euc_sim_h01  --recode --make-bed --pheno euc_sim_h01.pheno

./gemma-0.98.5-linux-static-AMD64 -bfile euc_sim_h10 -gk 1 -o euc_sim.gk -outdir .

./gemma-0.98.5-linux-static-AMD64 -bfile euc_sim_h10 -lmm 4 -k euc_sim.gk.cXX.txt -o euc_sim_h10.gemma -outdir .
./gemma-0.98.5-linux-static-AMD64 -bfile euc_sim_h01 -lmm 4 -k euc_sim.gk.cXX.txt -o euc_sim_h01.gemma -outdir .

./gemma-0.98.5-linux-static-AMD64 -bfile euc_sim_h10 -bslmm 1 -k euc_sim.gk.cXX.txt -o euc_sim_h10.bslmm -outdir .

