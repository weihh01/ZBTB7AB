cLoops2 pre -f DMSO_rep1.bedpe.gz -o DMSO_rep1 -p 40 
cLoops2 pre -f DMSO_rep2.bedpe.gz -o DMSO_rep2 -p 40 
cLoops2 pre -f DMSO_rep3.bedpe.gz -o DMSO_rep3 -p 40 

cLoops2 pre -f DMSO_rep1.bedpe.gz,DMSO_rep2.bedpe.gz,DMSO_rep3.bedpe.gz -o DMSO -p 40

cLoops2 pre -f dTAG_rep1.bedpe.gz -o dTAG_rep1 -p 40 
cLoops2 pre -f dTAG_rep2.bedpe.gz -o dTAG_rep2 -p 40 
cLoops2 pre -f dTAG_rep3.bedpe.gz -o dTAG_rep3 -p 40 

cLoops2 pre -f dTAG_rep1.bedpe.gz,dTAG_rep2.bedpe.gz,dTAG_rep3.bedpe.gz -o dTAG -p 40


cLoops2 estRes -d DMSO -o DMSO -bs 5000,10000 -p 40
cLoops2 estRes -d dTAG -o dTAG -bs 5000,10000 -p 40


cLoops2 estDis -d DMSO -o DMSO -bs 10000 -plot -p 40
cLoops2 estDis -d dTAG -o dTAG -bs 10000 -plot -p 40


cLoops2 callPeaks -d DMSO -o DMSO -eps 50,100 -minPts 10 -mcut 1000 -split -p 40
cLoops2 callPeaks -d dTAG -o dTAG -eps 50,100 -minPts 10 -mcut 1000 -split -p 40


cLoops2 callDiffLoops -tloop dTAG_loops.txt -td dTAG -cloop DMSO_loops.txt -cd DMSO -o dTAG_DMSO -j -w -p 40

cLoops2 agg -d dTAG -o 7A_dTAG_spe_dTAG -loops ../../CTCF/7A/cloops_rep1_rep2/dTAG_DMSO_dTAG_specific_dloops.txt -loop_norm -loop_vmin 0 -loop_vmax 1 -p 40
cLoops2 agg -d DMSO -o 7A_dTAG_spe_DMSO -loops ../../CTCF/7A/cloops_rep1_rep2/dTAG_DMSO_dTAG_specific_dloops.txt -loop_norm -loop_vmin 0 -loop_vmax 1 -p 40
