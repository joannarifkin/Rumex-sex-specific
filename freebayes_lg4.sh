#!/bin/bash

#freebayes version:  v1.1.0-60-gc15b070 

/ohta1/stephen.wright/freebayes/bin/freebayes \
  --fasta-reference /ohta1/bianca.sacchi/rhast_remap/genome/REF_LA.fa \
  --region LG4 \
	 /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_7.TM1_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_8.TM2_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_9.TM3_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_10.TM4_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_11.TM5_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_12.TM6_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_1.TF1_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_2.TF2_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_3.TF3_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_4.TF4_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_5.TF5_.bam /ohta1/bianca.sacchi/rhast_remap/JoshPop/STAR/Picard/AnalysisReady_6.TF6_.bam \
  --vcf JoshPopTX_LG4.vcf   \
~
