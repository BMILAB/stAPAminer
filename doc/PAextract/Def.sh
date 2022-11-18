#input_bam=$1
rdef_file=$1
rscAPA_file=$2
sample_name=$3
bam=$4
echo $sample_name
#保留unique map的结果
if [ ! -d "./temp" ]; then
  mkdir ./temp
fi
if [ ! -d "./log" ]; then
  mkdir ./log
fi
if [ ! -d "./DefSAF" ]; then
  mkdir ./DefSAF
fi
if [ ! -d "./NARSAF" ]; then
  mkdir ./NARSAF
fi
if [ ! -d "./scPA" ]; then
  mkdir ./scPA
fi
if [ ! -d "./DefExp" ]; then
  mkdir ./DefExp
fi
echo "保留unique map"
samtools view -@ 48 -h -F 256 -bS ${bam} > ./temp/filterBam.bam
samtools index -@ 24  ./temp/filterBam.bam
echo "根据位置去重"
# 10x
#umi_tools dedup -I ./temp/filterBam.bam -S ./temp/dedup_${sample_name}.bam \
#        --method=unique --extract-umi-method=tag \
#         --umi-tag=UB --cell-tag=CB

# no 10x
umi_tools dedup -I ./temp/filterBam.bam -S ./temp/dedup_${sample_name}.bam \ --method=unique
#切分BAM文件
echo "按方向切分BAM"
samtools view -@ 48 -h -F 0x10 -bS ./temp/dedup_${sample_name}.bam > ./temp/dedup.forward.bam
samtools view -@ 48 -h -f 0x10 -bS ./temp/dedup_${sample_name}.bam > ./temp/dedup.reverse.bam
samtools sort -@ 48 -o ./temp/dedup_${sample_name}.forward.sorted.bam ./temp/dedup.forward.bam && samtools index -@ 48 ./temp/dedup_${sample_name}.forward.sorted.bam
samtools sort -@ 48 -o ./temp/dedup_${sample_name}.reverse.sorted.bam ./temp/dedup.reverse.bam && samtools index -@ 48 ./temp/dedup_${sample_name}.reverse.sorted.bam 
echo "开始计算peak"
Rscript ${rdef_file} ./ temp dedup_${sample_name}.reverse.sorted.bam dedup_${sample_name}.forward.sorted.bam ./DefSAF/Def_${sample_name}.saf
echo "去重 定量"
featureCounts -a ./DefSAF/Def_${sample_name}.saf -F SAF -t exon -s 1 -M --largestOverlap -o ./temp/${sample_name}_peak_assigned -R BAM -T 24 ./temp/dedup_${sample_name}.bam > ./log/${sample_name}.FC.logfiles
samtools sort -@ 24 ./temp/dedup_${sample_name}.bam.featureCounts.bam -o ./temp/assigned_sorted.bam && samtools index -@ 24 ./temp/assigned_sorted.bam
umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --method=unique --per-cell -I ./temp/assigned_sorted.bam -S ./DefExp/${sample_name}.counts.tsv.gz > ./log/${sample_name}.count.logfiles
echo "计算scAPA位点"
Rscript ${rscAPA_file} ./temp/dedup_${sample_name}.bam ./scPA ${sample_name}
