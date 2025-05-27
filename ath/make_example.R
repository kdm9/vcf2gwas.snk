library(tidyverse)


ft16 = read_csv("https://arapheno.1001genomes.org/rest/phenotype/30/values.csv")

acc_1k1g = request("https://tools.1001genomes.org//vcfsubset/data/1135.json") %>%
    req_perform() %>%
    resp_body_json()

acclist = as.integer(unlist(do.call("rbind", acc_1k1g$rows)[,1]))

pheno = ft16 %>%
    transmute(indiv=accession_id, ft16=phenotype_value) %>%
    filter(!is.na(ft16), indiv%in% acclist) %>%
    write_tsv("ft16.csv", na="") %>%
    glimpse()

library(httr2)

r =  request('https://tools.1001genomes.org/api/v1/vcfsubset/') |>
    req_headers(
        Origin="https://tools.1001genomes.org",
        Referer="https://tools.1001genomes.org/vcfsubset/",
        Host="tools.1001genomes.org",
        Accept="*/*"
    ) |>
    req_body_form(
        strains=paste(sprintf("%d", pheno$indiv), collapse=","),
        regions="5:3000000-4000000",
        type="fullgenome",
        format="vcf"
    ) |>
    req_progress() |>
    req_perform("ath")

system("bcftools view -i 'MAC>5 & F_MISSING < 0.2' ath.vcf --write-index -Oz -o ath_filt-MAC5-MISS20.vcf.gz")

