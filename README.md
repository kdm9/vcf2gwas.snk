# vcf2gwas.snk

A partial snakemake reimplementation of Vogt's vcf2gwas

# "documentation"

I am working on more complete docs, but for now please refer to the following.


The core functionality of this is parallelised GEMMA across phenotypes, including all file format conversion hocus-pocus. 

To run it, you need four things:

- The snakefile: `vcf2gwas.snk`. This is not user-serviceable unless you want to extend what vcf2gwas does
- A config file, e.g. `config.yml`. This provides all user settings to vcf2gwas. The example is based on *A. thaliana* flowering time. This can be called anything, just make sure you give the correct file to `--configfile` in the snakemake command below. See the comments in `config.yml` for further explanation of config options.
- A vcf. It should be standards compliant and indexed. If in doubt, run it through `bcftools view --write-index` to make sure.
- A phenotype file. This should be a csv or tsv, whose first column is called `indiv`, and is the individual name. Each phenotype should be a separate column, with a column name in the header. Try to use names that are valid R/python column/variable names, i.e. no spaces and avoid punctuation.

Sample names must match exactly between the phenotype file and the vcf. If the names don't match, you must manually rename them by editing the tsv or using `bcftools annotate` to change the names in the VCF.

To run vcf2gwas, edit (a copy of) the config file to your preferred settings. Then, run `vcf2gwas` like so:

    # ONCE, install the binaries using:
    bash install_binaries.sh

    # to run the pipeline:
    snakemake --snakefile vcf2gwas.snk --cores all --configfile config.yml

