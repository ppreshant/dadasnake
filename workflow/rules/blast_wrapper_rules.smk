# Blast related rules from the snakemake wrapper. Reference :
# https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast.html?highlight=blast

# name of the custom database file # put into config?
DB_SOURCE_FILENAME = '9-species-16S' # skip the .fasta extension


# run blastn : nucleotides ; adapted to run inside Dadasnake workflow -- directories and such
# Source : https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/blast/blastn.html
rule blast_nucleotide:
    input:
        query = "sequenceTables/all.seqs.fasta",  # general command: {sample}.fasta
        blastdb=multiext("../DBs/blast/" + DB_SOURCE_FILENAME + '/database',
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
        format="6 qseqid sseqid evalue",
        extra=""
    wrapper:
        "v1.21.0/bio/blast/blastn"


# Make a custom BLAST database from fasta file

rule blast_makedatabase_nucleotide:
    input:
        fasta="../reference_sequences/" + DB_SOURCE_FILENAME + '.fasta' # path relative to dadasnake or the output directory
    output:
        DB_OUTPUTS = multiext("../DBs/blast/" + DB_SOURCE_FILENAME + '/database',
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
            )
    log:
        "../DBs/blast/" + DB_SOURCE_FILENAME + "/" + DB_SOURCE_FILENAME + ".log"
    params:
        "-input_type fasta -blastdb_version 5 -parse_seqids"
    wrapper:
        "v1.21.0/bio/blast/makeblastdb"
# run independently inside DBs folder with
# snakemake -s ../dadasnake/workflow/rules/blast_wrapper_rules.smk -c 4 --use-conda blast_makedatabase_nucleotide

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
