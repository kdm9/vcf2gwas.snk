#!/bin/bash
set -eu
mkdir -p .opt/bin .opt/tmp
cd .opt/tmp
curl -LSsO https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip
unzip -q plink_linux_x86_64_20241022.zip
mv plink  ../bin

curl -LSsO https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
gunzip -f gemma-0.98.5-linux-static-AMD64.gz
mv gemma-0.98.5-linux-static-AMD64 ../bin/gemma
chmod +x ../bin/gemma

curl -LSsO https://github.com/precimed/simu/releases/download/v0.9.2/simu_linux
mv simu_linux ../bin/simu
chmod +x ../bin/simu

curl -LSsO https://yanglab.westlake.edu.cn/software/gcta/bin/gcta-1.94.3-linux-kernel-3-x86_64.zip
unzip -q gcta-1.94.3-linux-kernel-3-x86_64.zip
mv  gcta-1.94.3-linux-kernel-3-x86_64/gcta64 ../bin/gcta64

curl -LSsO https://github.com/shenwei356/csvtk/releases/download/v0.32.0/csvtk_linux_amd64.tar.gz
tar xf csvtk_linux_amd64.tar.gz
mv csvtk ../bin/

cd ..
rm -rf tmp
ls bin
