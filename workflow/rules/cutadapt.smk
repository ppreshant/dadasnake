def get_fastq(wildcards):
    fastqs = config['raw_directory'] + "/" + samples.loc[(wildcards.library,wildcards.run), ["r1_file", "r2_file"]].dropna()
    return fastqs

def get_lib_perRunAndSample(wildcards,prefix,suffix):
    return prefix+samples.loc[(samples['run']==wildcards.run) & (samples['sample']==wildcards.sample), "library"].unique()+suffix


localrules: primers_control

rule primers_control: # confirms that the rules relating to primers are done
    input:
        expand("preprocessing/{samples.run}/{samples.sample}.{direction}.fastq.gz", samples=samples.itertuples(), direction=["fwd","rvs"]),
        "reporting/readNumbers.tsv",
        "reporting/primerNumbers_perSample.tsv"
    output:
        "primers.done"
    shell:
        """
        touch {output}
        """

rule combine_or_rename:
    input:
        "reporting/primerNumbers_perLibrary.tsv",
        files = lambda wildcards: get_lib_perRunAndSample(wildcards,"preprocessing/{run}/",".{direction}.fastq.gz")
    output:
        "preprocessing/{run}/{sample}.{direction}.fastq.gz"
    wildcard_constraints:
        direction="(fwd|rvs)",
        sample='|'.join(samples['sample'])
    threads: 1
    log: "logs/combine_or_rename.{run}.{sample}.{direction}.log"
    resources:
        runtime="01:00:00",
        mem=config['normalMem']
    run:
        if len(input) > 2:
            shell("cat {input.files} > {output}")
        else:
            shell("mv {input.files} {output}")

# Count the raw reads
rule input_numbers:
    input:
        "reporting/sample_table.tsv",
        expand("{raw_directory}/{file}", file=samples.r1_file,raw_directory=RAW),
        expand("{raw_directory}/{file}", file=samples.r2_file,raw_directory=RAW)
    output:
        report("reporting/readNumbers.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "raw",
        raw_directory = RAW
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    conda: ENVDIR + "dada2_env.yml"
    log: "logs/countInputReads.log"
    script:
        SCRIPTSDIR+"report_readNumbers.R"


# Count reads after primer processing
rule primer_numbers:
    input:
        "reporting/readNumbers.tsv",
        expand("preprocessing/{samples.run}/{samples.library}.{direction}.fastq.gz", samples=samples.itertuples(), direction=["fwd","rvs"])
    output:
        report("reporting/primerNumbers_perLibrary.tsv",category="Reads"),
        report("reporting/primerNumbers_perSample.tsv",category="Reads")
    threads: 1
    params:
        currentStep = "primers"
    resources:
        runtime="12:00:00",
        mem=config['normalMem']
    log: "logs/countPrimerReads.log"
    conda: ENVDIR + "dada2_env.yml"
    script:
        SCRIPTSDIR+"report_readNumbers.R"


# run cutadapt to trim primers
if config['sequencing_direction'] == "fwd_1":
    rule cut_primer_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fwd.fastq.gz",
            "preprocessing/{run}/{library}.rvs.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming forward primer is in read 1. {config[primers][fwd][sequence]}"
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} 'XXXXXX')
            FWD_RC=`echo {config[primers][fwd][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            RVS_RC=`echo {config[primers][rvs][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`

            cutadapt -g {config[primers][fwd][sequence]} -G {config[primers][rvs][sequence]} \
            {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
            -O {config[primer_cutting][overlap]} \
            -m 1:1 --pair-filter={config[primer_cutting][filter_if_not_match]} \
            -j {threads} -e {config[primer_cutting][perc_mismatch]} --trimmed-only \
            -o $TMPD/{wildcards.library}.fwd.fastq.gz -p $TMPD/{wildcards.library}.rvs.fastq.gz {input} &> {log}

            cutadapt -a $RVS_RC -A $FWD_RC \
             {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
             -m 1:1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             -o {output[0]} -p {output[1]} $TMPD/{wildcards.library}.fwd.fastq.gz $TMPD/{wildcards.library}.rvs.fastq.gz >> {log} 2>&1
             """

elif config['sequencing_direction'] == "rvs_1":
    rule cut_primers_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fwd.fastq.gz",
            "preprocessing/{run}/{library}.rvs.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Assuming forward primer is in read 2."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} "XXXXXX")

            FWD_RC=`echo {config[primers][fwd][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            RVS_RC=`echo {config[primers][rvs][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`

            cutadapt -g {config[primers][rvs][sequence]} -G {config[primers][fwd][sequence]} \
             {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
              -O {config[primer_cutting][overlap]} \
              -m 1:1 --pair-filter={config[primer_cutting][filter_if_not_match]} \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} --trimmed-only \
             -o $TMPD/{wildcards.library}.rvs.fastq.gz -p $TMPD/{wildcards.library}.fwd.fastq.gz {input} &> {log}

            cutadapt -a $RVS_RC -A $FWD_RC \
             {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
              -m 1:1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             -o {output[0]} -p {output[1]} $TMPD/{wildcards.library}.fwd.fastq.gz $TMPD/{wildcards.library}.rvs.fastq.gz >> {log}  2>&1
            """

else:
    rule cut_primers_both:
        input:
            get_fastq
        output:
            "preprocessing/{run}/{library}.fwd.fastq.gz",
            "preprocessing/{run}/{library}.rvs.fastq.gz"
        threads: 1
        resources:
            runtime="12:00:00",
            mem=config['normalMem']
        conda: ENVDIR + "dadasnake_env.yml"
        log: "logs/cutadapt.{run}.{library}.log"
        message: "Running cutadapt on {input}. Searching for both  primers in both reads."
        shell:
            """
            TMPD=$(mktemp -d -t --tmpdir={TMPDIR} "XXXXXX")
            FWD_RC=`echo {config[primers][fwd][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`
            RVS_RC=`echo {config[primers][rvs][sequence]} | tr '[ATUGCYRSWKMBDHNatugcyrswkbdhvn]' '[TAACGRYSWMKVHDBNtaacgryswmkvhdbn]' |rev`

            cutadapt -g {config[primers][fwd][sequence]} -G {config[primers][rvs][sequence]} \
             {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
              -O {config[primer_cutting][overlap]} \
              -m 1:1 --pair-filter={config[primer_cutting][filter_if_not_match]} \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             --untrimmed-output=$TMPD/{wildcards.library}.fwd_unt.fastq.gz --untrimmed-paired-output=$TMPD/{wildcards.library}.rvs_unt.fastq.gz \
             -o $TMPD/{wildcards.library}.fwd.fastq.gz -p $TMPD/{wildcards.library}.rvs.fastq.gz {input} &> {log}

            cutadapt -a $RVS_RC -A $FWD_RC \
             {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
              -m 1:1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             -o {output[0]} -p {output[1]} $TMPD/{wildcards.library}.fwd.fastq.gz $TMPD/{wildcards.library}.rvs.fastq.gz >> {log} 2>&1

            cutadapt -g {config[primers][rvs][sequence]} -G {config[primers][fwd][sequence]} \
             {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
              -O {config[primer_cutting][overlap]} \
              -m 1:1 --pair-filter={config[primer_cutting][filter_if_not_match]} \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} --trimmed-only \
             -o $TMPD/{wildcards.library}.rvs_tr2.fastq.gz -p $TMPD/{wildcards.library}.fwd_tr2.fastq.gz \
             $TMPD/{wildcards.library}.fwd_unt.fastq.gz $TMPD/{wildcards.library}.rvs_unt.fastq.gz >> {log} 2>&1

            cutadapt -a $RVS_RC -A $FWD_RC \
             {config[primer_cutting][indels]} -n {config[primer_cutting][count]} \
              -m 1:1 \
             -j {threads} -e {config[primer_cutting][perc_mismatch]} \
             -o $TMPD/{wildcards.library}.fwd.final.fastq.gz -p $TMPD/{wildcards.library}.rvs.final.fastq.gz \
             $TMPD/{wildcards.library}.fwd_tr2.fastq.gz $TMPD/{wildcards.library}.rvs_tr2.fastq.gz >> {log} 2>&1

            cat $TMPD/{wildcards.library}.fwd.final.fastq.gz >> {output[0]}
            cat $TMPD/{wildcards.library}.rvs.final.fastq.gz >> {output[1]}
            """
