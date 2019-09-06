#! /bin/bash -l

#proj=facility
email=iryna.shutava@umu.se

#module load bioinfo-tools
#module load MEME

## default args
in=/mnt/picea/projects/arabidopsis/cbellini/ninja-motif-seach/MEME/DE_novel_motifs/SubGenome.fa
out=/mnt/picea/projects/arabidopsis/cbellini/ninja-motif-seach/MEME/DE_novel_motifs/Results
bkgrModel=/mnt/picea/projects/arabidopsis/cbellini/ninja-motif-seach/MEME/DE_novel_motifs/bckgrModel_Full
cpu=8

## create the out dir
#if [ ! -d $out ]; then
#    mkdir -p $out
#fi

## prepare
sbatch -n 8 --array=0-8 --mail-user=$email -e $out/SubGenome_DE_MEME.fa_%a.err \
-o $out/SubGenome_DE_MEME.fa_%a.out ../UPSCb-common/pipeline/runMeme.sh $in $out $bkgrModel

#sbatch -A facility -p all -n 8 --array=1-667 -c $cpu -e $out/SubGenome_DE_MEME.fa_%a.err -o $out/SubGenome_DE_MEME.fa_%a.out \
#../UPSCb-common/pipeline/runBlastPlus2.sh -f 5 -p $cpu blastx $in $inx $out $bkgrModel
