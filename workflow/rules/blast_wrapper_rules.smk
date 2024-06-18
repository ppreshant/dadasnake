# Blast related rules from the snakemake wrapper. Reference :
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast.html?highlight=blast

# name of the custom database file
DB_SOURCE_FILENAME = '6-species-S090' # skip the .fasta extension
# TODO: put into config if integrating this into the dadasnake pipeline


# run blastn : nucleotides ; adapted to run inside Dadasnake workflow -- directories and such
# Source : https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast/blastn.html
rule blast_nucleotide:
    input:
        query = "sequenceTables/all.seqs.fasta",  # general command: {sample}.fasta
        blastdb=multiext("../../DBs/blast/" + DB_SOURCE_FILENAME + '/database',
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
    output:
        "mapping/all.seqs.blast.txt"
    log:
        "logs/all.seqs.blast.log"
    threads:
        2
    params:
        # Usable options and specifiers for the different output formats are listed here:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast/blastn.html.
        format="6 qseqid sseqid evalue mismatch length",
        extra=""
    wrapper:
        "v1.21.0/bio/blast/blastn"
# note: currently one entry for each sequence is provided, this can be changed by max_target_seqs
# ref: https://www.ncbi.nlm.nih.gov/books/NBK279684/. I am using an R script to filter the top hit

# Simplify BLAST output. Select top hit per query : with custom R script
rule simplify_blast:
    input:
        "mapping/all.seqs.blast.txt",
        summary_table = "sequenceTables/all.seqTab.tsv"
    log:
        "logs/all.seqs.blast.log"
    threads:
        2
    output:
        "mapping/top_all.seqs.blast.tsv",
        summary_out = "sequenceTables/all.seqTab_custom-blast.tsv"
    conda:
        ENVDIR + "tidyverse_env.yml"
    script:
        SCRIPTSDIR + "simplify_blast.R"


# Collate BLAST output and plot. Merge ASV counts mapping to same organism
rule collate_plot_blast:
    input:
        "sequenceTables/all.seqTab_custom-blast.tsv"
    log:
        "logs/cleanup.blast.log"
    threads:
        1
    output:
        collated_summary = "sequenceTables/collated.seqTab_custom-blast.tsv",
        plot_linear = 'mapping/plot_collated_blast.pdf',
        plot_log = 'mapping/plot-log_collated_blast.pdf'
    conda:
        ENVDIR + "tidyverse_env.yml"
    script:
        SCRIPTSDIR + "collate_plot_blast.R"


# Make a custom BLAST database from fasta file
# Note: DB_SOURCE_FILENAME needs to be entered in the beginning of this file
# paths here are relative to the output directory of dadasnake workflow
rule blast_makedatabase_nucleotide:
    input:
        fasta="../../reference_sequences/" + DB_SOURCE_FILENAME + '.fasta' # path relative to dadasnake or the output directory
    output:
        DB_OUTPUTS = multiext("../../DBs/blast/" + DB_SOURCE_FILENAME + '/database',
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
            )
    log:
        "../../DBs/blast/" + DB_SOURCE_FILENAME + "/" + DB_SOURCE_FILENAME + ".log"
    params:
        "-input_type fasta -blastdb_version 5 -parse_seqids"
    wrapper:
        "v1.21.0/bio/blast/makeblastdb"
# run with this command, inside dadasnake/ in the snakemake_conda env
# snakemake --configfile config/config.RAM_16S.yaml --use-conda --conda-prefix $PWD/conda --cores 3 blast_makedatabase_nucleotide
# or below command runs it from working dir of dadasnake/dadasnake_output
# snakemake -c 3 --use-conda blast_makedatabase_nucleotide


# remove the specific blast database : folder and subcontents.. USED ONLY FOR TESTING PURPOSES
rule clear_database:
    params:  DB_SOURCE_FILENAME
    shell: 'rm -r ../DBs/blast/{DB_SOURCE_FILENAME}'


# ------------------

# protein database
rule blast_makedatabase_protein:
    input:
        fasta="protein/{protein}.fasta"
    output:
        multiext("../blast-results/{protein}.fasta",
            ".pdb",
            ".phr",
            ".pin",
            ".pot",
            ".psq",
            ".ptf",
            ".pto"
        )
    log:
        "logs/{protein}.log"
    params:
        "-input_type fasta -blastdb_version 5"
    wrapper:
        "v1.21.0/bio/blast/makeblastdb"


# run blastp for proteins
