mkdir .opt
mkdir .opt/bin .opt/tmp
cd .opt/tmp
curl -LSsO https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip
unzip plink_linux_x86_64_20241022.zip
mv plink  ../bin

curl -LSsO https://github.com/genetics-statistics/GEMMA/releases/download/v0.98.5/gemma-0.98.5-linux-static-AMD64.gz
gunzip -f gemma-0.98.5-linux-static-AMD64.gz
mv gemma-0.98.5-linux-static-AMD64 ../bin/gemma
chmod +x ../bin/gemma

cd ..
rm -rf tmp
ls bin
