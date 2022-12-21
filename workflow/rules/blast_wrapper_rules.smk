# Blast related rules from the snakemake wrapper. Reference :
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast.html?highlight=blast

rule blast_makedatabase_nucleotide:
    input:
        fasta="9-species-16S.fasta" # path relative to dadasnake or the output directory
    output:
        multiext("blast/9-species-16S.fasta",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
    log:
        "../DBs/blast/9-species-16S.log"
    params:
        "-input_type fasta -blastdb_version 5 -parse_seqids"
    wrapper:
        "v1.21.0/bio/blast/makeblastdb"
# run independantly inside DBs folder with
# snakemake -s ../dadasnake/workflow/rules/blast_wrapper_rules.smk -c 4 --use-conda blast_makedatabase_nucleotide


# run blastn : nucleotides
# Source : https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast/blastn.html
rule blast_nucleotide:
    input:
        query = "{sample}.fasta",
        blastdb=multiext("blastdb/blastdb",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
    output:
        "{sample}.blast.txt"
    log:
        "logs/{sample}.blast.log"
    threads:
        2
    params:
        # Usable options and specifiers for the different output formats are listed here:
        # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast/blastn.html.
        format="6 qseqid sseqid evalue",
        extra=""
    wrapper:
        "v1.21.0/bio/blast/blastn"



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
