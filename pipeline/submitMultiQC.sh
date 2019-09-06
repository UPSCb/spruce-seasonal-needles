#!/bin/bash -l

set -ex

# environment variables
proj=u2015030
mail="thri0013@gapps.umu.se"
in=/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030/fastqc
out=/mnt/picea/projects/spruce/sjansson/seasonal-needles/u2015030/multiqc

#load necessary modules
module load bioinfo-tools multiqc

#checks if UPSCb environment variable is defined in our git UPSCb checkout directory
if [ -z ../UPSCb-common ]; then 
    echo "Set up the UPSCb env. var. to your Git UPSCb checkout dir." 
    exit 1
fi

#checks if input file exists
if [ ! -d $in ]; then
	echo "The first argument needs to be an existing analysis directory."
fi

#checks if output directory exist, if not creates it
if [ ! -d $out ]; then
  mkdir $out
fi

#also verify the options
sbatch -p all -w aspseq -A $proj -e $out.err -o $out.out --mail-user $mail ../UPSCb-common/pipeline/runMultiQC.sh $in $out

# -A, --account : Charge resources used by this job to specified account.
# -e,  --error=<filename pattern>:
# Instruct Slurm to connect the batch script's standard error directly to the file name specified in the "filename pattern". 
# By default both standard output and standard error are directed to the same file. 
# For job arrays, the default file name is "slurm-%A_%a.out", "%A" is replaced by the job ID and "%a" with the array index.
# -o, --output=<filename pattern>:
# Instruct Slurm to connect the batch script's standard output directly to the file name specified in the "filename pattern".

# -p, --partition 'all' can acces every server (picea,watson,aspseq) (partitions can be: node, core, all)
# -w choose aspseq as working server
