#! /bin/bash
#SBATCH --time=0-08:00:00
#SBATCH --mem-per-cpu=16G
#SBATCH --ntasks=2

source $CONDA_ACTIVATE deeptools

workDir=/mnt/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq
runName=diff_abund_3_canonical_noRRnoSP_moreSeq
regionData="autosomalArms"
chipData="publicChIPseq"
outDir=${workDir}/${runName}/custom/exploreChIP
mkdir -p $outDir


bwlistfile=${workDir}/${runName}/custom/exploreChIP/publicChIPseq.tsv

sampleNames=( $(cut -f 1 $bwlistfile) )
echo ${sampleNames[@]}
fileNames=( $(cut -f 3 $bwlistfile) )
echo ${fileNames[@]}


## bigwig files
bwPaths=${fileNames[@]}
bwNames=${sampleNames[@]}

## bed files 
upReg=( `ls ${workDir}/${runName}/custom/exploreChIP/bed/*autosomalArmsUp*.bed` )
upRegNames=( `basename -a ${upReg[@]} | awk -F"__" '{print $1}'` )
downReg=( `ls ${workDir}/${runName}/custom/exploreChIP/bed/*autosomalArmsDown*.bed` )
downRegNames=( `basename -a ${downReg[@]} | awk -F"__" '{print $1}'` )


echo "number of tasks: " $SLURM_NTASKS

## compute matrix
redoMatrix=true
if $redoMatrix; then 
      blacklist=https://github.com/Boyle-Lab/Blacklist/raw/master/lists/ce11-blacklist.v2.bed.gz
      blacklistFile=`basename ${blacklist%.gz}`
      if [ ! -f "$blacklistFile" ]; then
            wget $blacklist
            gunzip `basename $blacklist`
      fi
 
      
      computeMatrix scale-regions -S ${bwPaths[@]} \
            -R ${upReg[@]} \
            --beforeRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --afterRegionStartLength 1000 \
            --skipZeros \
            -o ${outDir}/matrix_${regionData}_up_${chipData}.mat.gz \
	      --outFileNameMatrix ${outDir}/matrix_${regionData}_up_${chipData}.tab \
            --blackListFileName $blacklistFile \
            --numberOfProcessors $SLURM_NTASKS \
            --verbose
      
      computeMatrix scale-regions -S ${bwPaths[@]} \
            -R ${downReg[@]} \
            --beforeRegionStartLength 1000 \
            --regionBodyLength 2000 \
            --afterRegionStartLength 1000 \
            --skipZeros \
            -o ${outDir}/matrix_${regionData}_down_${chipData}.mat.gz \
	      --outFileNameMatrix ${outDir}/matrix_${regionData}_down_${chipData}.tab \
            --blackListFileName $blacklistFile \
            --numberOfProcessors $SLURM_NTASKS \
            --verbose
fi

# multiBigwigSummary BED-file --bwfiles ${bwPaths[@]} --BED $activeEnhancers_fountain \
#     --outRawCounts matrix_enhVhistoneModEncode_avrRaw.tab -o matrix_enhVhistoneModEncode_avrRaw.npz \
#     --numberOfProcessors 4 

# maxs=( `awk 'NR > 1 { for (i = 4; i <= NF; i++) { if ($i > max[i]) max[i] = $i+0 } } END { for (i = 4; i <= NF; i++) print max[i]/2+0 }' matrix_enhVhistoneModEncode_avrRaw.tab` )
#maxs=( 30 `printf ' 2 %.0s' {1..21}` )
#mins=( `printf ' -2 %.0s' {1..22}` )
maxs=( 60 300 1000 500 10 200 800 200 1000 )
mins=( `printf ' 0 %.0s' {1..9}` )


## make plots
plotHeatmap -m ${outDir}/matrix_${regionData}_up_${chipData}.mat.gz \
      -out ${outDir}/${regionData}_up_${chipData}_heatmap.pdf \
      --colorMap Blues  \
      --startLabel "TSS" --endLabel "TES" \
      -y "" -x "Distance" \
      --regionsLabel ${upRegNames[@]}  \
      --plotTitle "Autosomal arm up regulated genes" \
      --samplesLabel ${bwNames[@]} \
      --sortRegions no \
      --zMin ${mins[@]} \
      --zMax ${maxs[@]} \
      --yMin ${mins[@]} \
      --yMax ${maxs[@]} 
# --sortUsingSamples 1 2 3 \

## make plots
plotHeatmap -m ${outDir}/matrix_${regionData}_down_${chipData}.mat.gz \
      -out ${outDir}/${regionData}_down_${chipData}_heatmap.pdf \
      --colorMap Blues  \
      --startLabel "TSS" --endLabel "TES" \
      -y "" -x "Distance" \
      --regionsLabel ${downRegNames[@]} \
      --plotTitle "Autosomal arm down regulated genes" \
      --samplesLabel ${bwNames[@]} \
      --sortRegions no \
      --zMin ${mins[@]} \
      --zMax ${maxs[@]} \
      --yMin ${mins[@]} \
      --yMax ${maxs[@]} 
#--sortUsingSamples 1 2 3 \

# plotProfile -m ${outDir}/matrix_${regionData}_${chipData}.mat.gz  \
#               -out ${outDir}/${regionData}_${chipData}_profile.png \
#               --numPlotsPerRow 5 \
#               --regionsLabel "upReg" "downReg" \
#               --startLabel "prom" --endLabel "" \
#               --colors "cyan" "magenta" "brown" "grey" \
#               --yMax ${maxs[@]} \
#               --yMin ${mins[@]} \
#               --samplesLabel ${bwNames[@]} \
# 	      --outFileNameData ${outDir}/${regionData}_${chipData}_profile.tab \
#               --plotTitle "Promoters of COH-1 regulated genes"
