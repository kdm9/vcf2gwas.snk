#!/usr/bin/env Rscript

library(argparser, quietly=TRUE)
options(argparser.delim="|")

ap = arg_parser("Plot a GEMMA lmm result")

ap = add_argument(ap, "--assocs", nargs="+", help=".assoc.txt file from a GEMMA -lmm run")
ap = add_argument(ap, "--pval", help="Which p-value to use", default="p_wald")
ap = add_argument(ap, "--output", help="Output file")
ap = add_argument(ap, "--causal-dir", help="directory with simu causal files")

args = parse_args(ap)

if (is.na(args$output)) {
    args$output = paste0("sim_plot.png")
}

library(readr, quietly=TRUE)
library(dplyr, quietly=TRUE)
library(tidyr, quietly=TRUE)
library(ggplot2, quietly=TRUE)

options(readr.show_col_types=F)
gwas = data.frame(file=unlist(strsplit(args$assocs, "|", fixed=T))) %>%
    mutate(
        data=purrr::map(file, function(x) {
            pheno=sub(".*/([^/]+).assoc.txt", "\\1", x, perl=T)
            simcausal = sprintf("%s/%s.1.causals", args$causal_dir, pheno)
            cau = read_tsv(simcausal) %>%
                transmute(chr=CHR, ps=POS, causal=T)
            read_tsv(x) %>%
                mutate(pheno=pheno) %>%
                rename(pval={!!args$pval}) %>%
                left_join(cau, by=join_by(chr, ps)) %>%
                mutate(causal=ifelse(is.na(causal), F, causal)) %>%
                arrange(causal, chr, ps)
    })) %>%
    unnest(data)

p = ggplot(gwas, aes(x=ps, y=-log10(pval))) +
    geom_point(aes(colour=causal, alpha=causal, size=causal)) +
    scale_alpha_manual(values=c(0.1, 1)) +
    scale_size_manual(values=c(0.1, 1)) +
    scale_colour_manual(values=c("black", "red")) +
    facet_grid(pheno~chr, scales="free_x", space="free_x") +
    theme_bw()
ggsave(args$output, plot=p, dpi=600, width=8, height=4)

