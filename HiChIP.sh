for i in 7A_CTCF_DMSO 7A_CTCF_dTAG 7B_CTCF_DMSO 7B_CTCF_dTAG;do
        trim_galore -q 20 -j 7 --phred33 --stringency 3 --length 20 -e 0.1 --paired ${i}_R*fq.gz -o ./
        cutadapt -a file:/home/whh/private/linker_micro-c/linker/MicroC_S-linker.fa -A file:/home/whh/private/linker_micro-c/linker/MicroC_S-linker.fa -j 10 --minimum-length=10 -o trimLk/${i}_R1_val_1_trimlinker.fq.gz -p trimLk/${i}_R2_val_2_trimlinker.fq.gz ${i}_R1_val_1.fq.gz ${i}_R2_val_2.fq.gz > trimLk/${i}.txt
        cutadapt -a file:/home/whh/private/linker_micro-c/linker/MicroC_S-linker.fa -A file:/home/whh/private/linker_micro-c/linker/MicroC_S-linker.fa -j 10 --minimum-length=10 -o trimLk/${i}_R1_val_1_trimlinker2.fq.gz -p trimLk/${i}_R2_val_2_trimlinker2.fq.gz trimLk/${i}_R1_val_1_trimlinker.fq.gz trimLk/${i}_R2_val_2_trimlinker.fq.gz > trimLk/${i}2.txt
        cutadapt -a file:/home/whh/private/linker_micro-c/linker/MicroC_S-linker_AT.fa -A file:/home/whh/private/linker_micro-c/linker/MicroC_S-linker_AT.fa -j 10 --minimum-length=10 -o trimLk/${i}_R1_val_1_trimlinker3.fq.gz -p trimLk/${i}_R2_val_2_trimlinker3.fq.gz trimLk/${i}_R1_val_1_trimlinker2.fq.gz trimLk/${i}_R2_val_2_trimlinker2.fq.gz > trimLk/${i}3.txt
	mkdir data_${i}
	cd data_${i}
	mkdir ${i}
	cd ..
	cp trimLk/${i}_R1_val_1_trimlinker3.fq.gz data_${i}/${i}/${i}_trimlinker3_1.fq.gz
	cp trimLk/${i}_R2_val_2_trimlinker3.fq.gz data_${i}/${i}/${i}_trimlinker3_2.fq.gz
	/opt/HiC-Pro-3.0.0/bin/HiC-Pro -c config-hicpro.txt -i data_${i}/ -o results_${i} -s mapping -s quality_checks
	/opt/HiC-Pro-3.0.0/bin/HiC-Pro -c config-hicpro.txt -i results_${i}/bowtie_results/bwt2 -o results_${i}/ -s proc_hic -s merge_persample -s quality_checks
	cd CTCF_ChIP
	bamToBed -bedpe -i ../results_${i}/bowtie_results/bwt2/${i}/${i}_trimlinker3_hg38XX.bwt2pairs.bam | sort -k1,1 -k2,3n -k4,4 -k5,6n -S 80G | awk '{if($1==$4)print$0}' OFS='\t' | awk '{if($2>$5)print $4 "\t" $5 "\t" $6 "\t" $1 "\t" $2 "\t" $3 "\t" $7 "\t" $8 "\t" $10 "\t" $9;else print $0}' | sort -k1,6 -k8,8 -k9,10 -u -S 90G | awk '{if($1!~/chr[CLMT]/ && $4!~/chr[CLMT]/) print $0}' > ${i}_trimlinker4_hg38XX.bedpe
	awk '{print $1"\t"int(($2+$3)/2-50)"\t"int(($2+$3)/2+50)"\n"$4"\t"int(($5+$6)/2-50)"\t"int(($5+$6)/2+50)}' ${i}_trimlinker4_hg38XX.bedpe > ${i}.bed
	num=$(cat ${i}.bed | awk 'END{print FNR}')
	sort -k1,1 -k2,2n -k3,3n -S 90G ${i}.bed | genomeCoverageBed -bg -i - -g /ssd/genome/hg38_chromsize.txt | awk -va=$num '{printf "%s\t%.0f\t%.0f\t%.2f\n",$1,$2,$3,$4/a*1000000}' | awk '{$4/=1;print}' OFS='\t' > ${i}.bdg
	bedGraphToBigWig ${i}.bdg /ssd/genome/hg38_chromsize.txt ${i}.bw
	cd ..
done
