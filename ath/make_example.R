library(tidyverse)
library(glue)
library(httr2)


ft16 = read_csv("https://arapheno.1001genomes.org/rest/phenotype/30/values.csv")

acc_1k1g = request("https://tools.1001genomes.org/vcfsubset/data/1135.json") %>%
    req_perform() %>%
    resp_body_json()

acclist = as.integer(unlist(do.call("rbind", acc_1k1g$rows)[,1]))

pheno = ft16 %>%
    transmute(indiv=accession_id, ft16=phenotype_value) %>%
    filter(!is.na(ft16), indiv%in% acclist) %>%
    write_tsv("ft16.tsv", na="") %>%
    glimpse()


#r =  request('https://tools.1001genomes.org/api/v1/vcfsubset/') |>
#    req_headers(
#        Origin="https://tools.1001genomes.org",
#        Referer="https://tools.1001genomes.org/vcfsubset/",
#        Host="tools.1001genomes.org",
#        Accept="*/*"
#    ) |>
#    req_body_form(
#        strains=paste(sprintf("%d", pheno$indiv), collapse=","),
#        regions="5:2000000-4000000",
#        type="fullgenome",
#        format="vcf"
#    ) |>
#    req_progress() |>
#    req_perform("ath.vcf")
#
samplescsv =paste(sprintf("%d", pheno$indiv), collapse=",")
region="5:2000000-4000000"
url="https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz"
writeLines(sprintf("%d", pheno$indiv), "samplelist.txt")

bed=expand.grid(chr=1:5, start=2000000, region=5000000) %>%
    write_tsv("regions_2m-5m-allchr.bed", col_names=F)
