# include configuration file
configfile: srcdir("config/config.default.yaml")

SCRIPTSDIR = srcdir("workflow/scripts/")
ENVDIR = srcdir("workflow/envs/")
ROOTDIR = srcdir("/")

include:
    "workflow/rules/get_config.smk"

workdir:
    OUTPUTDIR

f = open('full.config.yaml', 'w+')
yaml.dump(config, f, allow_unicode=True,default_flow_style=False)

def getThreads(max):
    if workflow.cores:
        realThreads = max if max <= workflow.cores else workflow.cores
    elif workflow.nodes:
        realThreads = max if max <= workflow.nodes else workflow.nodes
    else:
        realThreads = max
    return realThreads


# processing paired end Reads
if config['paired']:
    if 'primers' in STEPS:
        include:
            "workflow/rules/cutadapt.smk"
    else:
        include:
            "workflow/rules/copying.smk"
    if 'dada' in STEPS:
        if not config['dada']['pool']:
            if not config['big_data']:
                include:
                    "workflow/rules/dada.paired.smk"
            else:
                include:
                    "workflow/rules/bigdada.paired.smk"
        else:
            if config['dada']['pool']=="within_run":
                include:
                    "workflow/rules/dada.paired.runpool.smk"
            else:
                include:
                    "workflow/rules/dada.paired.pool.smk"

# processing non-paired end reads
else:
    if 'primers' in STEPS:
        include:
            "workflow/rules/cutadapt.single.smk"
    else:
        include:
            "workflow/rules/copying.single.smk"
    if 'dada' in STEPS:
        if not config['dada']['pool']:
            if not config['big_data']:
                include:
                    "workflow/rules/dada.single.smk"
            else:
                include:
                    "workflow/rules/bigdada.single.smk"
        else:
            if config['dada']['pool']=="within_run":
                include:
                    "workflow/rules/dada.single.runpool.smk"
            else:
                include:
                    "workflow/rules/dada.single.pool.smk"

# common processing for paired and non-paired end reads
if 'dada' in STEPS:
    if config['big_data']:
        include:
            "workflow/rules/bigdada.common.smk"
    else:
        include:
            "workflow/rules/dada.common.smk"
if 'taxonomy' in STEPS:
    if config['big_data']:
        include:
            "workflow/rules/bigtaxonomy.smk"
    else:
        include:
            "workflow/rules/taxonomy.smk"
if 'postprocessing' in STEPS:
    if config['final_table_filtering']['do']:
        if config['big_data']:
            include:
                "workflow/rules/bigpost.filtering.smk"
        else:
            include:
                "workflow/rules/post.filtering.smk"
    else:
        if config['big_data']:
            include:
                "workflow/rules/bigpost.no_filtering.smk"
        else:
            include:
                "workflow/rules/post.no_filtering.smk"


inputs = []
if 'primers' in STEPS:
    inputs.append('primers.done')
if 'dada' in STEPS:
    inputs.append('dada.done')
if 'taxonomy' in STEPS:
    inputs.append('taxonomy.done')
if 'postprocessing' in STEPS:
    inputs.append('postprocessing.done')
if config['hand_off']['biom'] and ( 'dada' in STEPS or 'taxonomy' in STEPS or 'postprocessing' in STEPS):
    inputs.append('sequenceTables/all.seqTab.biom')

if EMAIL == "":
    onsuccess:
        shell("mkdir -p job.errs.outs &>> logs/cleanup.log; ( mv dadasnake* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *stdout job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *log job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *logfile job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log")
else:
    onsuccess:
        shell('mkdir -p job.errs.outs &>> logs/cleanup.log; ( mv dadasnake* job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *stderr job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *stdout job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *log job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; ( mv *logfile job.errs.outs || touch job.errs.outs ) &>> logs/cleanup.log; echo "$(date) {config[sessionName]}" | mail -s "dadasnake finished" {EMAIL} ')
    onerror:
        shell('echo "$(date) {config[sessionName]}" | mail -s "dadasnake exited with error" {EMAIL} ')
    onstart:
        shell('echo "$(date) {config[sessionName]}" | mail -s "dadasnake started" {EMAIL} ')

localrules: ALL
#, SamplesPrint

# master command
rule ALL:
    input:
        inputs
    output:
        touch('workflow.done')

rule SamplesPrint:
    input:
        sam_path
    output:
        "reporting/sample_table.tsv"
    threads: 1
    log: "logs/printsamples.log"
    resources:
        runtime="01:00:00",
        mem=config['normalMem']
    run:
        samples.to_csv(path_or_buf=output[0],sep="\t",index=False,index_label=False)

DIRECTION_LIST = ['fwd', 'rvs'] # to reuse
import glob

# Custom rules : Executed separately from the workflow

# Make a fasta file of all the OTUs in the final taxonomy table tsv
rule taxonomyTable_to_fasta: # run from the datasnake directory within conda snakemake_env
    input: "sequenceTables/all.seqTab.tax.tsv"
    output: "OTUsequences.fa"
    message: "Saved OTUs into a single fasta file"
    shell: # Remove header, grab first two columns from input, make fasta
        r"""
        touch {output}
        cat {input} | tail +2 | awk -F'\t' '{{print ">" $1 "\n" $2}}' > {output}
        """

# Custom rule to clean or rerun specific segments of the workflow
# This rule should not run in the regular pipeline since it has no outputs.. But this can only be run once
rule clean_filtered_outputs:
    input: # identify the outputs of all the filter commands to delete
        glob.glob("filtered/*"), # filtered output files
        expand('stats/fastqc_filtered_{dir}/', dir = DIRECTION_LIST),  # stats of filtered samples
        expand("stats/multiqc_filtered_{dir}_report.html", dir = DIRECTION_LIST),
        expand("stats/QC_filtered.{run}.{dir}.pdf", run=samples.run.unique(), dir = DIRECTION_LIST), # PDF outputs
        glob.glob("reporting/filteredNumbers_per*"),  # reports
        "sequenceTables/all.seqTab.RDS", # inputs to data_control -- makes dada.done final
        "sequenceTables/all.seqs.fasta",
        "reporting/finalNumbers_perSample.tsv"
        #'dada.done' # the marker for the dada2 pipeline; useful to run subsequent steps ; run it alone
    message: "Successfully removed files {input}. \n You can rerun the filtering step now using these as targets"
    shell: "rm -r {input}"

# after clearing these intermediate files, run snakemake with these files as target:
# stats/multiqc_filtered_{fwd,rvs}_report.html stats/QC_filtered.1.{fwd,rvs}.pdf reporting/filteredNumbers_perLibrary.tsv

# include rules for BLAST wrapper
include:
    "workflow/rules/blast_wrapper_rules.smk"
