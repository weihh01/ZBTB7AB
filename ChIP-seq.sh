for i in ChIP-DMSO ChIP-Flavo;do
	trim_galore -q 20 -j 7 --phred33 --stringency 3 --length 20 -e 0.1 --paired ${i}_R*fq.gz -o ./
	bowtie2 -p 60 --no-unal --no-mixed  --no-discordant --dovetail --very-sensitive --score-min L,0,-0.4 -X 1000 -x /mnt/disk4/bowtie2/hg38XX -1 ${i}_R1_val_1.fq.gz -2 ${i}_R2_val_2.fq.gz | samtools view -bh -q 10 > ${i}_q10.bam
	samtools sort -@ 40 -o ${i}_q10.sort.bam ${i}_q10.bam
	java -jar /mnt/disk1/6/share/software/picard.jar MarkDuplicates -REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT -I ${i}_q10.sort.bam -O ${i}_q10_markdup_DIS10k.bam -M ${i}_DIS10k.metrics --OPTICAL_DUPLICATE_PIXEL_DISTANCE 10000 --MAX_OPTICAL_DUPLICATE_SET_SIZE 300000 --TAGGING_POLICY All > ${i}_DIS10k.log 2>&1
	samtools index ${i}_q10_markdup_DIS10k.bam
	bamCoverage --binSize 10 --normalizeUsing CPM --effectiveGenomeSize 2913022398 --minMappingQuality 30 --ignoreDuplicates --centerReads -b ${i}_q10_markdup_DIS10k.bam -o ${i}_CPM.bw -p 40
done
