from debian:12

run export DEB_FRONTEND=noninteractve && apt update && apt install -y python3-pip python3-venv python3-dev build-essential curl unzip r-cran-tidyverse && python3 -m pip install snakemake --break-system-packages && mkdir /wd && Rscript -e 'install.packages("argparser")'

add vcf2gwas.snk ath.sh install_binaries.sh plot_gemma.R /wd/
workdir /wd/
run bash -xe install_binaries.sh
add https://pub.gekkonid.com/2025-02_1k1g-subset/1k1g_chr1:1000000-2000000_miss10.vcf.gz.csi /wd/
add https://pub.gekkonid.com/2025-02_1k1g-subset/1k1g_chr1:1000000-2000000_miss10.vcf.gz /wd/
