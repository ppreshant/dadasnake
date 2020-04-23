#! /bin/bash -i

CONFIGFILE=$1
VARCONFIG=$2
DIR=$(dirname $VARCONFIG)
JNAME=$3

while read var val; do unset $var ; declare $var="$val" ; done < $VARCONFIG
if [ "$SNAKEMAKE_VIA_CONDA" = true ]; then
   CONDA_START="conda activate $DIR/snakemake_env"
   CONDA_END="conda deactivate"
else
   CONDA_START=""
   CONDA_END=""
fi

eval $LOADING_MODULES
eval $CONDA_START

snakemake -j 50 -s $DIR/Snakefile --cluster-config $DIR/dada_scripts/slurm.config.yaml --cluster "{cluster.call} {cluster.runtime}{params.runtime} {cluster.mem_per_cpu}{params.mem} {cluster.threads}{threads} {cluster.partition}" --configfile $CONFIGFILE --config sessionName=$JNAME --use-conda --conda-prefix $DIR/dada_env_common >> $JNAME.stdout 2>> $JNAME.stderr

snakemake -j 1 -s $DIR/Snakefile --report report.html --configfile $CONFIGFILE

eval $CONDA_END
