################ 0. Import packages ################


#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("ff")
#install.packages("ggplot2")
#install.packages("stringr")
#install.packages("plotly")
#install.packages("gridExtra")
#install.packages("fuzzyjoin")
#install.packages("car")
#install.packages("lme4")
#install.packages("glmmTMB")
#install.packages("DHARMa")
#install.packages("corrplot")    
#install.packages("bbmle") ## for AICtab
#install.packages("car")
#install.packages("pscl")
#install.packages("fitdistrplus")
#install.packages("MASS")
#install.packages("tweedie")
#install.packages("colorBlindness")
#install.packages("Hmisc")
#install.packages("tibble")
#install.packages("egg")
#install.packages("gtable")
#install.packages("ppcor")

library(dplyr)
library(tidyr)
library(ff)
library(ggplot2)
library(stringr)
library(plotly)
library(gridExtra)
library(lme4)
library(glmmTMB)
library(DHARMa)
library(corrplot)    
library("bbmle") ## for AICtab
source("http://highstat.com/Books/BGS/GAMM/RCodeP2/HighstatLibV6.R")
library(car)
library(pscl)
library(fitdistrplus)
library(MASS)
library(tweedie)
library(colorBlindness)
library(purrr)
library(Hmisc)
library(tibble)
library(grid)
library(egg)
library(gtable)
library(ppcor)
#library(fuzzyjoin)
options(scipen=999)
options(tibble.width = 10, tibble.print_max = 10, tibble.print_min = 10)

#Analyses in R 4.1.0

#If running on ohta use #ohta commands 

################ 1. Import and combine data ################

  ################ 1.1 Import and combine linkage maps and calculate crossover count ################
    
    #Import the final physically anchored linkage maps and combine them into a single data table

    
    # setwd("/ohta1/joanna.rifkin/Lep-Map/Rumex/Sex_specific_map") #ohta


    chr1<-read.table("Joanna_final/phys/order1.mapped", stringsAsFactors = F)
    colnames(chr1)<-c("LG", "bp_position", "cm_male","cm_female")
    chr2<-read.table("Joanna_final/phys/order2.mapped", stringsAsFactors = F)
    colnames(chr2)<-c("LG", "bp_position", "cm_male","cm_female")
    chr3<-read.table("Joanna_final/phys/order3.mapped", stringsAsFactors = F)
    colnames(chr3)<-c("LG", "bp_position", "cm_male","cm_female")
    chr4<-read.table("Joanna_final/phys/order4.mapped", stringsAsFactors = F)
    colnames(chr4)<-c("LG", "bp_position", "cm_male","cm_female")
    chr5<-read.table("Joanna_final/phys/order5.mapped", stringsAsFactors = F)
    colnames(chr5)<-c("LG", "bp_position", "cm_male","cm_female")

    #setwd("/ohta1/joanna.rifkin/HiCSNPData/Sex_specific_map/Windowed_analyses") #ohta

    #Use inverse Haldane mapping function to add crossovers (x 186 individuals) to LG maps
        
    chr1 <- chr1 %>%  mutate (., diff_male=c(0, diff(chr1$cm_male))) %>% mutate (., crossovers_male= round(186*0.5*(1 - exp(-diff_male/50)))          )  %>% mutate (., diff_female=c(0, diff(chr1$cm_female))) %>% mutate (., crossovers_female= round(186*0.5*(1 - exp(-diff_female/50))))

    chr2 <- chr2 %>%  mutate (., diff_male=c(0, diff(chr2$cm_male))) %>% mutate (., crossovers_male= round(186*0.5*(1 - exp(-diff_male/50)))          )  %>% mutate (., diff_female=c(0, diff(chr2$cm_female))) %>% mutate (., crossovers_female= round(186*0.5*(1 - exp(-diff_female/50))))
    
    chr3 <- chr3 %>%  mutate (., diff_male=c(0, diff(chr3$cm_male))) %>% mutate (., crossovers_male= round(186*0.5*(1 - exp(-diff_male/50)))          )  %>% mutate (., diff_female=c(0, diff(chr3$cm_female))) %>% mutate (., crossovers_female= round(186*0.5*(1 - exp(-diff_female/50))))
    
    chr4 <- chr4 %>%  mutate (., diff_male=c(0, diff(chr4$cm_male))) %>% mutate (., crossovers_male= round(186*0.5*(1 - exp(-diff_male/50)))          )  %>% mutate (., diff_female=c(0, diff(chr4$cm_female))) %>% mutate (., crossovers_female= round(186*0.5*(1 - exp(-diff_female/50))))

    chr5 <- chr5 %>%  mutate (., diff_male=c(0, diff(chr5$cm_male))) %>% mutate (., crossovers_male= round(186*0.5*(1 - exp(-diff_male/50)))          )  %>% mutate (., diff_female=c(0, diff(chr5$cm_female))) %>% mutate (., crossovers_female= round(186*0.5*(1 - exp(-diff_female/50))))
    
    
    
    fullmap<-bind_rows(chr1, chr2, chr3, chr4, chr5)

    rm(chr1, chr2, chr3, chr4, chr5)
    
    fullmap<-fullmap %>% rowwise %>% mutate(.,cm_sex_averaged=mean(c(cm_male, cm_female)))
    fullmap<-fullmap %>% rowwise %>% mutate(.,crossovers_sex_averaged=mean(c(crossovers_male, crossovers_female)))
    
    #write.csv(fullmap, "Windowed_analyses/fullmap_1.1_8-25.csv", quote = F, row.names = F)
    #write.csv(fullmap, "fullmap_1.1.csv", quote = F, row.names = F) #ohta
    
    #View(fullmap)
    
    
 
  ################ 1.2 Attach scaffold positions ################

    #Add SNP names from old assembly - for conversion between map versions
    # setwd("/ohta1/joanna.rifkin/Lep-Map/Rumex/Sex_specific_map") #ohta
    
    snplist<-read.table("Joanna_final/snps.txt", header = T, stringsAsFactors = F)
    snplist$index<-seq.int(nrow(snplist))
    snplist$sex_linked<-as.numeric(grepl("\\*", snplist$POS))
    snplist$scaff_pos_num<-as.numeric(gsub("*", "",  x=snplist$POS, fixed=TRUE))
    colnames(snplist)[1:2]<-c("scaffold","scaff_pos")
    snplist_to_LGs<-read.table("Joanna_final/snps_c.liftover2.lg", header = T, stringsAsFactors = F)
    snplist_with_LGs_and_original_scaffolds<-left_join(snplist, snplist_to_LGs, by=c("index"="X0"))
    full_map_with_original_SNPs<-left_join(snplist_with_LGs_and_original_scaffolds, fullmap, by=c("CHR"="LG","POS"="bp_position"))
    full_map_with_original_SNPs_unique<-distinct(full_map_with_original_SNPs) #There are some duplicated rows in the original map file (I don't remember quite why but I know this is a known but harmless quirk for LepMap)
    rm(snplist, snplist_to_LGs, snplist_with_LGs_and_original_scaffolds, fullmap, full_map_with_original_SNPs)
    #setwd("/ohta1/joanna.rifkin/HiCSNPData/Sex_specific_map/Windowed_analyses") #ohta
    #write.csv(full_map_with_original_SNPs_unique, "fullmap_1.2.csv", quote = F, row.names = F) #ohta
    #write.csv(full_map_with_original_SNPs_unique, "Windowed_analyses/fullmap_1.2_8-23.csv", quote = F, row.names = F)    
    
  ################ 1.3 Attach segregation distortion data ################
    
    #Add segregation distortion LOD scores, calculated from X2.awk and X3.awk, see Joanna_final/collect_distortion_rates.txt
    #Distortion is calculated only at recombining markers - i.e., many NAs are expected

    sex_specific_distortion<-read.table("Segregation_distortion/6-9-2021_sex_specific_distortion_full_map.txt")
    colnames(sex_specific_distortion)<-c("marker_index","file_index","male_distortion","female_distortion")
    Mendelian_distortion<-read.table("Segregation_distortion/6-9-2021_Mendelian_distortion_full_map.txt")
    colnames(Mendelian_distortion)<-c("marker_index","file_index","Mendelian_distortion")
    distortion_all<-left_join(sex_specific_distortion, Mendelian_distortion, by="marker_index")
    full_map_with_original_SNPs_unique_and_distortion<-left_join(full_map_with_original_SNPs_unique, distortion_all, by=c("index"="marker_index") )
    full_map_with_original_SNPs_unique_and_distortion<- full_map_with_original_SNPs_unique_and_distortion %>% separate(., scaffold, into = c(NA, "scaffold_numeric",NA), sep="_|\\;", remove = F)
    
    
    View(full_map_with_original_SNPs_unique_and_distortion %>% group_by(CHR) %>% summarise(., sum(male_distortion>3.841, na.rm=T))) #Male distortion - number of significantly distorted sites on different LGs
    View(full_map_with_original_SNPs_unique_and_distortion %>% group_by(CHR) %>% summarise(., sum(female_distortion>3.841, na.rm=T)))      
    View(full_map_with_original_SNPs_unique_and_distortion %>% group_by(CHR) %>% summarise(., sum(female_distortion>6.635, na.rm=T)))     
    
    
    
    
    
    rm(distortion_all, full_map_with_original_SNPs_unique, Mendelian_distortion, sex_specific_distortion)
    #View(full_map_with_original_SNPs_unique_and_distortion)
    #write.csv(full_map_with_original_SNPs_unique_and_distortion, "fullmap_1.3.csv", quote = F, row.names = F) #ohta
    #write.csv(full_map_with_original_SNPs_unique_and_distortion, "Windowed_analyses/fullmap_1.3_8-23.csv", quote = F, row.names = F)
    
    
    
    
  ################ 1.4. Window map ################        

    
    #Average variables by 1mb windows across genome
    
    fullmap<-full_map_with_original_SNPs_unique_and_distortion
    rm(full_map_with_original_SNPs_unique_and_distortion)
    
    str(fullmap)
    colnames(fullmap)[7]<-"LG"
    colnames(fullmap)[8]<-"bp_position"
    
    
    #View(windowed_map)
    #View(fullmap)


#    test<-c(1.1, 1.1, NA, 1.1, NA)
    
  #  length(unique(test[!is.na(test)]))-1
    
    #This windows the map
    windowed_map<-fullmap %>% 
      group_by(LG) %>% #first, group by LG
      mutate(position_window=bp_position%/%as.numeric(1000000)) %>% 
      #Then, assign each SNP a window (i.e., which 1mb chunk of the chromosome is it in)
      group_by(LG,position_window) %>%  
      #group by position windows
      add_tally() %>% 
   #   na.omit() %>%
      #count SNPs per window as we go by
      select_if(., is.numeric) %>% #select numeric data
      group_by(LG,position_window) %>% #group again
      #next, calculate amount of centimorgans in that window for males and females ...
     # mutate(., m_crossovers=length(unique(cm_male))-1)  %>%
      #mutate(., f_crossovers=length(unique(cm_female))-1)  %>%
      mutate(., cm_span_male=(max(cm_male, na.rm=T)-min(cm_male,na.rm=T)))%>%
      mutate(., cm_span_female=(max(cm_female, na.rm=T)-min(cm_female, na.rm=T)))%>%
      #followed by cm/mbp (can be adjusted for different sized windows)
      mutate(., cm_mb_male=(max(cm_male, na.rm=T)-min(cm_male, na.rm=T))/(1000000)) %>% #I think there's actually an extra /1000000 in here, troubeshoot if using cm/mb
      mutate(., cm_mb_female=(max(cm_female, na.rm=T)-min(cm_female, na.rm=T))/(1000000)) %>% 
      mutate(., cm_mb_sex_averaged=(max(cm_sex_averaged, na.rm=T)-min(cm_sex_averaged, na.rm=T))/(1000000)) %>% 
      #add the cm-mb correlation 
      mutate(., windowed_correlation=cor(cm_male,cm_female, use="pairwise.complete.obs")) %>% 
      mutate_if(is.numeric, list(~na_if(., Inf))) %>% 
      #We're feeding SNP positions without CM positions into max, which results in negative infinities, so just reset all those to NAs
      summarise_at(vars(n,cm_male,cm_female,cm_sex_averaged,cm_span_male,cm_span_female,cm_mb_male,cm_mb_female,cm_mb_sex_averaged,crossovers_male,crossovers_female,crossovers_sex_averaged,bp_position, windowed_correlation, male_distortion, female_distortion, Mendelian_distortion),c(mean, sum), na.rm=T) %>% #calculate means for columns of interest
      mutate(., window_end = (position_window+1)*1000000 ) #and translate position windows back into more human readable numbers
    #   View(windowed_map)
    
    warnings()
    #Fix remaining NaNs and infinities
    
    windowed_map[windowed_map == -Inf] <- NA
    windowed_map[windowed_map == Inf] <- NA
    
    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    
    windowed_map[is.nan(windowed_map)] <- NA

    #Select and rename column names
    colnames(windowed_map) <- c("LG", "position_window", "n_mean", "cm_male_mean", "cm_female_mean", "cm_sex_averaged_mean", "cm_span_male_mean", "cm_span_female_mean", "cm_mb_male_mean", "cm_mb_female_mean", "cm_mb_sex_averaged_mean", "crossovers_male_mean", "crossovers_female_mean", "crossovers_sex_averaged_mean", "bp_position_mean", "windowed_correlation_mean", "male_distortion_mean", "female_distortion_mean", "Mendelian_distortion_mean", "n_sum", "cm_male_sum", "cm_female_sum", "cm_sex_averaged_sum", "cm_span_male_sum", "cm_span_female_sum", "cm_mb_male_sum", "cm_mb_female_sum", "cm_mb_sex_averaged_sum", "crossovers_male_sum", "crossovers_female_sum", "crossovers_sex_averaged_sum", "bp_position_sum", "windowed_correlation_sum", "male_distortion_sum", "female_distortion_sum", "Mendelian_distortion_sum", "window_end") 
    
 #   test<- windowed_map %>% select(., c("LG", "position_window", "cm_span_male_sum","cm_male_sum"))
    #    View(test)

    #View(windowed_map)    
    warnings()
    
    #write.csv(fullmap, "Windowed_analyses/fullmap_1.4_9-10-2021.csv", quote = F, row.names = F)
    #write.csv(windowed_map, "Windowed_analyses/windowed_map_1.4_8-30-2021.csv", quote = F, row.names = F)
        
    rm(fullmap)
################ 2. Window and attach gene data ################    

    ################ 2.1 Attach gene start and end positions ################
    
    #Add positions and identities of annotated genes
    

    
        ################ 2.1.1 Yunchen's new annotation ################
    
        #Import new annotation, filtered to include only gene heads (rather than the full list of mrnas or exons) 
        #Collate with version of annotation filtered to exclude anything with an exon that overlaps with a TE from Solomiya's TE annotation

        yunchen_annotation<-read.table("Windowed_analyses/gene_heads_only_REF_LA_genes.gff", stringsAsFactors = F, sep="\t", quote="")
        #Add quote character because of this weird issue: https://kbroman.org/blog/2017/08/08/eof-within-quoted-string/
        

        colnames(yunchen_annotation)<-c("yunchen_LG","yunchen_source","yunchen_type","yunchen_start","yunchen_end","NA","yunchen_strand","NA","yunchen_ID")
        
        yunchen_annotation<-yunchen_annotation[,c(1:5,7,9)]
        
        #Add gene lengths to the annotation
        yunchen_annotation$yunchen_length=abs(yunchen_annotation$yunchen_end-yunchen_annotation$yunchen_start)
        
        #Import exons that do not intersect TEs
        
        yunchen_annotation_filtered_exons<-read.table("Windowed_analyses/TE_filtered_REF_LA_exons.gff", stringsAsFactors = F, sep="\t")
      #  View(yunchen_annotation_exons)
        
        colnames(yunchen_annotation_filtered_exons)<-c("yunchen_LG","yunchen_source","yunchen_type","yunchen_start","yunchen_end","NA","yunchen_strand","NA","yunchen_ID")
        
        yunchen_annotation_filtered_exons<-yunchen_annotation_filtered_exons[,c(1:5,7,9)]
        yunchen_annotation_filtered_exons<-yunchen_annotation_filtered_exons  %>% separate(., yunchen_ID, into = c(NA,"yunchen_ID",NA,NA,NA,NA,NA,NA), sep = "\\;|\\=|\\:|\\-") 
 #       length(unique(yunchen_annotation$yunchen_ID))
  #      length(unique(yunchen_annotation_filtered_exons$yunchen_ID))
       # write.table(unique(yunchen_annotation_filtered_exons$yunchen_ID), "filtered_exon_list.txt", quote = F, row.names = F)

        #View(yunchen_annotation)
        
        #remove anything annotated we didn't already get that's clearly some kind of transposon 
        #Starting with 59121 genes in the annotation
        yunchen_strict_genes <- yunchen_annotation[!str_detect(yunchen_annotation$yunchen_ID, "[Tt]ranspos*"), ]  #49509
        #yunchen_strict_genes <- yunchen_strict_genes[!str_detect(yunchen_strict_genes$yunchen_ID, "unknown"), ] #Not necessary if TE filtered via TE annotation
        yunchen_strict_genes <- yunchen_strict_genes[!str_detect(yunchen_strict_genes$yunchen_ID, "[Rr]ibonuclease H"), ] #46626
        yunchen_strict_genes <- yunchen_strict_genes[!str_detect(yunchen_strict_genes$yunchen_ID, "[Pp]ol poly"), ] #46626 - total overlap with ribonucleaseH
        yunchen_strict_genes <- yunchen_strict_genes[!str_detect(yunchen_strict_genes$yunchen_ID, "[Mm]itochondri"), ] #44030 
        yunchen_strict_genes <- yunchen_strict_genes[!str_detect(yunchen_strict_genes$yunchen_ID, "[Cc]hloroplast"), ] #42388
        yunchen_strict_genes <- yunchen_strict_genes[!str_detect(yunchen_strict_genes$yunchen_ID, "[Rr]etrovirus"), ] #42388
        
      #  yunchen_strict_genes<-yunchen_strict_genes[,c(1,4:5,10)]
        yunchen_strict_genes<-yunchen_strict_genes  %>% separate(., yunchen_ID, into = c(NA,"yunchen_ID",NA,NA,NA,NA,NA,NA), sep = "\\;|\\=") 
        #Separate ID to match other IDs AFTER filtering by gene description
        
        yunchen_strict_genes<-subset(yunchen_strict_genes, yunchen_strict_genes$yunchen_ID %in% yunchen_annotation_filtered_exons$yunchen_ID) #30641 genes do not overlap TEs and meet our criteria
        

        
        ### Save position index of filtered gene set for later
        yunchen_position_index<-yunchen_strict_genes[,c(1,4:5,7:8)]
        

        #Window annotation for joining to map
                        
        windowed_yunchen_strict_genes<-yunchen_strict_genes %>% 
          group_by(yunchen_LG) %>%
          mutate(position_window=yunchen_start%/%as.numeric(1000000)) %>% 
          group_by(yunchen_LG,position_window) %>% 
          add_tally()%>% 
          select_if(., is.numeric) %>%
          group_by(yunchen_LG,position_window) %>%
          summarise_at(vars(n,yunchen_length), .funs = c(mean,sum), na.rm=T)
        
        #View(windowed_yunchen_strict_genes)
        
        windowed_yunchen_strict_genes<-windowed_yunchen_strict_genes[,c(1:3,6)]
        colnames(windowed_yunchen_strict_genes)<-  c("Yunchen_LG_1", "position_window", "n_genes_Yunchen","sum_Length_Yunchen")
        
        windowed_map_with_busco_and_yunchen_annotation<-left_join(windowed_map_with_busco, windowed_yunchen_strict_genes,by=c("LG"="Yunchen_LG_1", "position_window"="position_window"))
        #View(windowed_map_with_busco_and_yunchen_annotation)
        #View(windowed_yunchen_strict_genes)
        
        windowed_map_with_busco_and_yunchen_annotation$n_genes_Yunchen[is.na(windowed_map_with_busco_and_yunchen_annotation$n_genes_Yunchen)] <- 0 #Convert all NAs for number of buscos in windows to 0s because windows with no BUSCOs were not included in the GFF 
        
        windowed_map_with_busco_and_yunchen_annotation$sum_Length_Yunchen[is.na(windowed_map_with_busco_and_yunchen_annotation$sum_Length_Yunchen)] <- 0 #Convert all NAs for busco length sum in windows to 0s
        
                
        #write.csv(windowed_map_with_busco_and_yunchen_annotation, "Windowed_analyses/fullmap_2.1.3_8-23-2021.csv", quote = F, row.names = F)

        
        ## Make filtered gff for upload to COGE ##
        
        raw_annotation<-read.table("Windowed_analyses/TE_filtered_exon_filtered_REF_LA_genes.gff", stringsAsFactors = F, sep="\t", quote="")
        
        raw_annotation <- raw_annotation %>% separate(., V9, into = c(NA,"yunchen_ID",NA,NA,NA,NA,NA,NA), sep = "\\;|\\=|\\:|\\-", remove=F)
        
        filtered_annotation<-raw_annotation[raw_annotation$yunchen_ID %in%  yunchen_strict_genes$yunchen_ID,]
        
        filtered_annotation<-dplyr::select(filtered_annotation, -yunchen_ID)
        
        write.table(filtered_annotation, "Rumex_hastatulus_filtered_annotation.gff", row.names = F, col.names = F, quote = F, eol="\n", sep="\t")
        
        
        
        rm(windowed_yunchen_strict_genes)
        rm(yunchen_annotation)
        rm(yunchen_annotation_filtered_exons)
        rm(yunchen_strict_genes)
        rm(windowed_map_with_busco)
        rm(filtered_annotation)
        rm(raw_annotation)
        
        
 
    
    ################ 2.1.2 Attach TE start and end positions ################

    #Add positions and identities of annotated TEs    

        # /ohta1/joanna.rifkin/Lep-Map/Rumex/Sex_specific_map/Joanna_final/Windowed_analyses
       
        #First section run on HPC cluster for memory reasons
        
        TE_no_overlaps<-read.table("/ohta1/solomiya.hnatovska/genome2/EDTA/rumexEDTArun/shortenedheaders.REF_LA.fa.mod.EDTA.anno/shortenedheaders.REF_LA.fa.mod.EDTA.TEanno.split.gff3", stringsAsFactors = F, sep="\t")
        colnames(TE_no_overlaps)<-c("LG", "source", "sequence_ontology", "start", "end", "score", "strand", "phase", "attributes")
        unique(TE_no_overlaps$sequence_ontology)
        
        #[1] "Gypsy_LTR_retrotransposon"    "repeat_region"
        #[3] "LTR_retrotransposon"          "Copia_LTR_retrotransposon"
        #[5] "Mutator_TIR_transposon"       "target_site_duplication"
        #[7] "long_terminal_repeat"         "helitron"
        #[9] "Tc1_Mariner_TIR_transposon"   "PIF_Harbinger_TIR_transposon"
        #[11] "LINE_element"                 "CACTA_TIR_transposon"
        #[13] "hAT_TIR_transposon"
        
        TE_no_overlaps$TE_overlap_length<-abs(TE_no_overlaps$end-TE_no_overlaps$start)
        
        
        #Window, grouping by TE type as well as LG and position
        
        windowed_TE_no_overlaps<-TE_no_overlaps %>% 
          group_by(LG) %>%
          mutate(position_window=start%/%as.numeric(1000000)) %>% 
          group_by(LG,position_window,sequence_ontology) %>% 
          add_tally()%>% 
          select_if(., is.numeric) %>%
          group_by(LG,position_window,sequence_ontology) %>%
          summarise_at(vars(n,TE_overlap_length), .funs = c(mean,sum), na.rm=T)
        
        
        
        windowed_TE_no_overlaps<-windowed_TE_no_overlaps[,c(1:4,7:8)]
        
        colnames(windowed_TE_no_overlaps)<-c("LG", "position_window","sequence_ontology","N_TEs","sum_length_TEs","window_end")
        
        reshaped_windowed_TE_no_overlaps <- windowed_TE_no_overlaps %>% pivot_wider(names_from = sequence_ontology, values_from = c(N_TEs, sum_length_TEs))
        
        #write.csv(reshaped_windowed_TE_no_overlaps, "windowed_TE_no_overlaps_by_type.csv", quote = F, row.names = F)
        
        
        
        
        # Desktop part#
        #Open on desktop
        reshaped_windowed_TE_no_overlaps<-read.csv("Windowed_analyses/windowed_TE_no_overlaps_by_type.csv")
        
        reshaped_windowed_TE_no_overlaps[is.na(reshaped_windowed_TE_no_overlaps)] <- 0 #Convert all NAs in windows to 0s
        
        colnames(reshaped_windowed_TE_no_overlaps)
        
        reshaped_windowed_TE_no_overlaps$N_TEs_all<-rowSums(reshaped_windowed_TE_no_overlaps[,c(4:14)], na.rm = T)
        
        reshaped_windowed_TE_no_overlaps$sum_length_TEs_all<-rowSums(reshaped_windowed_TE_no_overlaps[,c(15:25)], na.rm = T)
        
        #Summarize by TE type
        
        #DNA transposons / Class 2
        
        #Tc1_Mariner_TIR_transposon
        #CACTA_TIR_transposon
        #hAT_TIR_transposon
        #PIF_Harbinger_TIR_transposon
        #Mutator_TIR_transposon
        #helitron
        
        reshaped_windowed_TE_no_overlaps$N_DNA_TEs_all<-rowSums(reshaped_windowed_TE_no_overlaps[,c(4,7,8,11,12,14)], na.rm = T)
        
        reshaped_windowed_TE_no_overlaps$sum_length_DNA_TEs_all<-rowSums(reshaped_windowed_TE_no_overlaps[,c(15,18,19,22,23,25)], na.rm = T)
        
        
        #RNA transposons / Class 1
        
        #N_TEs_Copia_LTR_retrotransposon
        #N_TEs_Gypsy_LTR_retrotransposon
        #N_TEs_LTR_retrotransposon
        #N_TEs_LINE_element
        
        
        reshaped_windowed_TE_no_overlaps$N_RNA_TEs_all<-rowSums(reshaped_windowed_TE_no_overlaps[,c(5,6,9,10)], na.rm = T)
        
        reshaped_windowed_TE_no_overlaps$sum_length_RNA_TEs_all<-rowSums(reshaped_windowed_TE_no_overlaps[,c(16,17,20,21)], na.rm = T)
        
        #View(reshaped_windowed_TE_no_overlaps)


        
        windowed_map_with_annotations_and_TEs<-left_join(windowed_map_with_busco_and_yunchen_annotation, reshaped_windowed_TE_no_overlaps,by=c("LG"="LG", "position_window"="position_window"))
        #View(windowed_map_with_annotations_and_TEs)
        
        
        #Set NAs to 0 because in this case absence of evidence IS evidence of absence
      
        windowed_map_with_annotations_and_TEs <- windowed_map_with_annotations_and_TEs %>% 
          mutate(
            N_TEs_CACTA_TIR_transposon = ifelse(is.na(N_TEs_CACTA_TIR_transposon), 0, N_TEs_CACTA_TIR_transposon),
            N_TEs_Copia_LTR_retrotransposon = ifelse(is.na(N_TEs_Copia_LTR_retrotransposon), 0, N_TEs_Copia_LTR_retrotransposon),
            N_TEs_Gypsy_LTR_retrotransposon = ifelse(is.na(N_TEs_Gypsy_LTR_retrotransposon), 0, N_TEs_Gypsy_LTR_retrotransposon),
            N_TEs_hAT_TIR_transposon = ifelse(is.na(N_TEs_hAT_TIR_transposon), 0, N_TEs_hAT_TIR_transposon),
            N_TEs_helitron = ifelse(is.na(N_TEs_helitron), 0, N_TEs_helitron),
            N_TEs_LINE_element = ifelse(is.na(N_TEs_LINE_element), 0, N_TEs_LINE_element),
            N_TEs_LTR_retrotransposon = ifelse(is.na(N_TEs_LTR_retrotransposon), 0, N_TEs_LTR_retrotransposon),
            N_TEs_Mutator_TIR_transposon = ifelse(is.na(N_TEs_Mutator_TIR_transposon), 0, N_TEs_Mutator_TIR_transposon),
            N_TEs_PIF_Harbinger_TIR_transposon = ifelse(is.na(N_TEs_PIF_Harbinger_TIR_transposon), 0, N_TEs_PIF_Harbinger_TIR_transposon),
            N_TEs_repeat_region = ifelse(is.na(N_TEs_repeat_region), 0, N_TEs_repeat_region),
            N_TEs_Tc1_Mariner_TIR_transposon = ifelse(is.na(N_TEs_Tc1_Mariner_TIR_transposon), 0, N_TEs_Tc1_Mariner_TIR_transposon),
            sum_length_TEs_CACTA_TIR_transposon = ifelse(is.na(sum_length_TEs_CACTA_TIR_transposon), 0, sum_length_TEs_CACTA_TIR_transposon),
            sum_length_TEs_Copia_LTR_retrotransposon = ifelse(is.na(sum_length_TEs_Copia_LTR_retrotransposon), 0, sum_length_TEs_Copia_LTR_retrotransposon),
            sum_length_TEs_Gypsy_LTR_retrotransposon = ifelse(is.na(sum_length_TEs_Gypsy_LTR_retrotransposon), 0, sum_length_TEs_Gypsy_LTR_retrotransposon),
            sum_length_TEs_hAT_TIR_transposon = ifelse(is.na(sum_length_TEs_hAT_TIR_transposon), 0, sum_length_TEs_hAT_TIR_transposon),
            sum_length_TEs_helitron = ifelse(is.na(sum_length_TEs_helitron), 0, sum_length_TEs_helitron),
            sum_length_TEs_LINE_element = ifelse(is.na(sum_length_TEs_LINE_element), 0, sum_length_TEs_LINE_element),
            sum_length_TEs_LTR_retrotransposon = ifelse(is.na(sum_length_TEs_LTR_retrotransposon), 0, sum_length_TEs_LTR_retrotransposon),
            sum_length_TEs_Mutator_TIR_transposon = ifelse(is.na(sum_length_TEs_Mutator_TIR_transposon), 0, sum_length_TEs_Mutator_TIR_transposon),
            sum_length_TEs_PIF_Harbinger_TIR_transposon = ifelse(is.na(sum_length_TEs_PIF_Harbinger_TIR_transposon), 0, sum_length_TEs_PIF_Harbinger_TIR_transposon),
            sum_length_TEs_repeat_region = ifelse(is.na(sum_length_TEs_repeat_region), 0, sum_length_TEs_repeat_region),
            sum_length_TEs_Tc1_Mariner_TIR_transposon = ifelse(is.na(sum_length_TEs_Tc1_Mariner_TIR_transposon), 0, sum_length_TEs_Tc1_Mariner_TIR_transposon),
            N_TEs_all = ifelse(is.na(N_TEs_all), 0, N_TEs_all),
            sum_length_TEs_all = ifelse(is.na(sum_length_TEs_all), 0, sum_length_TEs_all),
            N_DNA_TEs_all = ifelse(is.na(N_DNA_TEs_all), 0, N_DNA_TEs_all),
            sum_length_DNA_TEs_all = ifelse(is.na(sum_length_DNA_TEs_all), 0, sum_length_DNA_TEs_all),
            N_RNA_TEs_all = ifelse(is.na(N_RNA_TEs_all), 0, N_RNA_TEs_all),
            sum_length_RNA_TEs_all = ifelse(is.na(sum_length_RNA_TEs_all), 0, sum_length_RNA_TEs_all) 
            )

        colnames(windowed_map_with_annotations_and_TEs)
        #View(windowed_map_with_annotations_and_TEs)
        rm(reshaped_windowed_TE_no_overlaps)
        rm(windowed_map_with_busco_and_yunchen_annotation)
        
 #       plot(windowed_map_with_annotations_and_TEs$cm_mb_female, windowed_map_with_annotations_and_TEs$sum_length_TEs_all)

 #       plot(windowed_map_with_annotations_and_TEs$male_distortion, windowed_map_with_annotations_and_TEs$sum_length_TEs_all)
        
        
    ################ 3. Window and attach expression and DE status ################
    
    #Add results of differential expression analysis to each section of genome    
    getwd()
    
    #windowed_map_with_annotations_and_TEs<-read.csv("Windowed_analyses/fullmap_2.1.4_8-11.csv")
    
    
    ################ 3.1 Number of expressed genes ################
    
    #Use position index generated from Yunchen's annotation above    
        

   # rm(yunchen_annotation)
    
    ################ 3.1.1 Number of expressed genes in leaf ################
    
    expressed_leaf<-read.csv("Windowed_analyses/expressed_genes_leaf.csv")
    
    colnames(expressed_leaf)<-c("Index", "yunchen_ID")
    
    expressed_leaf$expressed_leaf<-1
    
    expressed_leaf<-expressed_leaf[,c(2:3)]
    
    expressed_leaf<-expressed_leaf  %>% separate(., yunchen_ID, into = c("yunchen_ID",NA), sep = "\\-") 
    
    expression<-left_join(yunchen_position_index, expressed_leaf, by="yunchen_ID")
    
    expression$expressed_leaf[is.na(expression$expressed_leaf)] <- 0 

    expression
    
   # View(expression)
    rm(expressed_leaf)
    rm(yunchen_position_index)
    
    ################ 3.1.2 Number of expressed genes in leaf ################
    
    expressed_flower<-read.csv("Windowed_analyses/expressed_genes_flower.csv")
    
    colnames(expressed_flower)<-c("Index", "yunchen_ID")
    
    expressed_flower$expressed_flower<-1
    
    expressed_flower<-expressed_flower[,c(2:3)]
    
    expressed_flower<-expressed_flower  %>% separate(., yunchen_ID, into = c("yunchen_ID",NA), sep = "\\-") 
    
    expression<-left_join(expression, expressed_flower, by="yunchen_ID")

    expression$expressed_flower[is.na(expression$expressed_flower)] <- 0 
    
    
    rm(expressed_flower)
    
    ################ 3.2 Number of differentially expressed genes ################
    
    ################ 3.2.1 Genes that are sex-biased in leaves ################
    
    sex_biased_leaf<-read.csv("Windowed_analyses/sex_biased_leaf_TX.csv")

    sex_biased_leaf<-subset(sex_biased_leaf, sex_biased_leaf$padj<=0.01)
    
    colnames(sex_biased_leaf)[4]<-"leaf_sex_log2FoldChange"

    sex_biased_leaf<-sex_biased_leaf %>%
      mutate(
        leaf_directional_bias = case_when(
          Bias == "MaleBias" ~ 1,
          Bias == "FemaleBias" ~ -1,
        )
      )    
    
    sex_biased_leaf$leaf_male_bias<-as.numeric(sex_biased_leaf$Bias=="MaleBias")
    
    sex_biased_leaf$leaf_female_bias<-as.numeric(sex_biased_leaf$Bias=="FemaleBias")
    
    sex_biased_leaf<-sex_biased_leaf  %>% separate(., GeneID, into = c("yunchen_ID",NA), sep = "\\-") 
    
    colnames(sex_biased_leaf)
    
    sex_biased_leaf<-sex_biased_leaf[,c(2,4,13:15)]

    expression<-left_join(expression, sex_biased_leaf, by="yunchen_ID")
    
    expression$leaf_male_bias[is.na(expression$leaf_male_bias)] <- 0 
    expression$leaf_female_bias[is.na(expression$leaf_female_bias)] <- 0 

    
    
    
   # View(sex_biased_leaf)
    rm(sex_biased_leaf)
    
    ################ 3.2.2 Genes that are sex-biased in flowerbuds ################

    
    sex_biased_flowerbud<-read.csv("Windowed_analyses/sex_biased_flowerbud_TX.csv")
    
    colnames(sex_biased_flowerbud)[4]<-"flowerbud_sex_log2FoldChange"
    
    sex_biased_flowerbud<-subset(sex_biased_flowerbud, sex_biased_flowerbud$padj<=0.01)
    
    sex_biased_flowerbud<-sex_biased_flowerbud %>%
      mutate(
        flower_directional_bias = case_when(
          Bias == "MaleBias" ~ 1,
          Bias == "FemaleBias" ~ -1,
        )
      )    
    
    sex_biased_flowerbud$flower_male_bias<-as.numeric(sex_biased_flowerbud$Bias=="MaleBias")
    
    sex_biased_flowerbud$flower_female_bias<-as.numeric(sex_biased_flowerbud$Bias=="FemaleBias")
    
    sex_biased_flowerbud<-sex_biased_flowerbud  %>% separate(., GeneID, into = c("yunchen_ID",NA), sep = "\\-") 
    
    colnames(sex_biased_flowerbud)
    
    sex_biased_flowerbud<-sex_biased_flowerbud[,c(2,4,13:15)]
    
    expression<-left_join(expression, sex_biased_flowerbud, by="yunchen_ID")

    
    expression$flower_male_bias[is.na(expression$flower_male_bias)] <- 0 
    expression$flower_female_bias[is.na(expression$flower_female_bias)] <- 0 
    
    
        
    rm(sex_biased_flowerbud)

        
    ################ 3.2.3 Genes that are pollen-biased ################

    pollen_biased<-read.csv("Windowed_analyses/pollen_biased_TX.csv")

    colnames(pollen_biased)[3]<-"pollen_log2FoldChange"
    
    pollen_biased<-subset(pollen_biased, pollen_biased$padj<=0.01)
    
    pollen_biased$pollen_biased<-1
    
    colnames(pollen_biased)
    
    pollen_biased<-pollen_biased  %>% separate(., gene, into = c("yunchen_ID",NA), sep = "\\-") 
    
    pollen_biased<-pollen_biased[,c(1,3,8)]
    
    expression<-left_join(expression, pollen_biased, by="yunchen_ID")

    expression$pollen_biased[is.na(expression$pollen_biased)] <- 0 
    
            
    rm(pollen_biased)

    
    ################ 3.2.4 Genes that are pollen tube-biased ################

    pollen_tube_biased<-read.csv("Windowed_analyses/pollenTube_biased_TX.csv")
    
    colnames(pollen_tube_biased)[3]<-"pollen_tube_log2FoldChange"
    
    pollen_tube_biased<-subset(pollen_tube_biased, pollen_tube_biased$padj<=0.01)
    
    pollen_tube_biased$pollen_tube_biased<-1
    
    pollen_tube_biased<-pollen_tube_biased  %>% separate(., gene, into = c("yunchen_ID",NA), sep = "\\-") 
    
    pollen_tube_biased<-pollen_tube_biased[,c(1,3,8)]
    
    expression<-left_join(expression, pollen_tube_biased, by="yunchen_ID")

    
    expression$pollen_tube_biased[is.na(expression$pollen_tube_biased)] <- 0 
    
        
    rm(pollen_tube_biased)
    
    
    
    ################ 3.2.5 Summarise expression bias ################
    
    View(expression %>% group_by(yunchen_LG) %>% count(leaf_male_bias, leaf_female_bias))
    View(expression %>% group_by(yunchen_LG) %>% count(flower_male_bias, flower_female_bias))

    View(expression %>% group_by(yunchen_LG) %>% count(pollen_biased, pollen_tube_biased))
    
        
    
        ################ 3.3 Window expression ################

    
    
    
    
        
    expression$genes_length<-abs(expression$yunchen_end-expression$yunchen_start)
    
    colnames(expression)
    
    windowed_expression<-expression %>% 
      group_by(yunchen_LG) %>%
      mutate(position_window=yunchen_start%/%as.numeric(1000000)) %>% 
      group_by(yunchen_LG,position_window) %>% 
      add_tally()%>% 
      select_if(., is.numeric) %>%
      group_by(yunchen_LG,position_window) %>%
      summarise_at(vars(
        expressed_leaf,expressed_flower,
        leaf_sex_log2FoldChange, leaf_directional_bias, leaf_male_bias,leaf_female_bias,
        flowerbud_sex_log2FoldChange, flower_directional_bias, flower_male_bias, flower_female_bias,
        pollen_log2FoldChange, pollen_biased,
        pollen_tube_log2FoldChange, pollen_tube_biased,
        genes_length
        ),
                   .funs = c(mean,sum), na.rm=T) %>%
      mutate_if(is.numeric, list(~na_if(., Inf))) %>%
      mutate_if(is.numeric, list(~na_if(., -Inf)))

   # View(windowed_expression)

    colnames(windowed_expression)<-c("yunchen_LG", "position_window", "expressed_leaf_mean", "expressed_flower_mean", "leaf_sex_log2FoldChange_mean", "leaf_directional_bias_mean", "leaf_male_bias_mean", "leaf_female_bias_mean", "flowerbud_sex_log2FoldChange_mean", "flower_directional_bias_mean", "flower_male_bias_mean", "flower_female_bias_mean", "pollen_log2FoldChange_mean", "pollen_biased_mean", "pollen_tube_log2FoldChange_mean", "pollen_tube_biased_mean", "genes_length_mean", "expressed_leaf_sum", "expressed_flower_sum", "leaf_sex_log2FoldChange_sum", "leaf_directional_bias_sum", "leaf_male_bias_sum", "leaf_female_bias_sum", "flowerbud_sex_log2FoldChange_sum", "flower_directional_bias_sum", "flower_male_bias_sum", "flower_female_bias_sum", "pollen_log2FoldChange_sum", "pollen_biased_sum","pollen_tube_log2FoldChange_sum", "pollen_tube_biased_sum","genes_length_sum")
    
    windowed_expression<-windowed_expression %>% 
      dplyr::select(., c(
        yunchen_LG, position_window,  #positional info
        expressed_leaf_sum, expressed_flower_sum, #number of expressed genes
        leaf_sex_log2FoldChange_mean, leaf_directional_bias_mean, leaf_male_bias_sum, leaf_female_bias_sum, #leaf sex biased expression
        flowerbud_sex_log2FoldChange_mean, flower_directional_bias_mean,flower_male_bias_sum, flower_female_bias_sum, #flowerbud sex biased expression
        pollen_log2FoldChange_mean, pollen_biased_sum, #pollen biased expression
        pollen_tube_log2FoldChange_mean, pollen_tube_biased_sum #pollen tube biased expression 
      ))

    is.nan.data.frame <- function(x)
      do.call(cbind, lapply(x, is.nan))
    windowed_expression[is.nan(windowed_expression)] <- NA
    
        

    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    
    length(unique(windowed_expression$yunchen_LG))
    
    length(unique(windowed_map_with_annotations_and_TEs$LG))
    
    setdiff(unique(windowed_expression$yunchen_LG),unique(windowed_map_with_annotations_and_TEs$LG))

 #   test<-subset(windowed_expression, windowed_expression$yunchen_LG %in% LGs)
    
   # test2<-subset(windowed_map_with_annotations_and_TEs, windowed_map_with_annotations_and_TEs$LG %in% LGs)

    complete_windowed_dataset<-left_join(windowed_map_with_annotations_and_TEs, windowed_expression, by=c("LG"="yunchen_LG", "position_window"))
    
    View(complete_windowed_dataset)
    
    rm(expression)
    rm(windowed_expression)
    
    
    #complete_windowed_dataset[,c(16:17,19:20,23:46, 48:49,52:53,56:57,59,61)]
     
    
    
    #write.csv(complete_windowed_dataset, "Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv",quote = F, row.names = F)
    
    
    
    
################ 3.4 Summaries of raw and windowed data ################
    
################ 3.4.1 Consecutive non-recombining windows ################
    #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
    
    #use the "rle" function to identify the lengths of consecutive windows with no crossovers
View(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG3"))
    
    #LG1
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG1")$crossovers_male_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG1")$crossovers_male_mean))))
    
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG1")$crossovers_female_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG1")$crossovers_female_mean))))
    
    
    #LG2
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG2")$crossovers_male_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG2")$crossovers_male_mean))))
    
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG2")$crossovers_female_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG2")$crossovers_female_mean))))
    
    
    #LG3
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG3")$crossovers_male_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG3")$crossovers_male_mean))))
    
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG3")$crossovers_female_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG3")$crossovers_female_mean))))
    
    
    #LG4
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG4")$crossovers_male_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG4")$crossovers_male_mean))))
    
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG4")$crossovers_female_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG4")$crossovers_female_mean))))
    
    
    #LG5
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG5")$crossovers_male_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG5")$crossovers_male_mean))))
    
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG5")$crossovers_female_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG5")$crossovers_female_mean))))
    
    
    
    
    complete_windowed_dataset$cm_mb_male_very_low<-(complete_windowed_dataset$cm_mb_male_mean<0.01)
    
    
    #LG1
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG1")$cm_mb_male_very_low))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG1")$crossovers_male_mean))))
    
    View(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG1")$crossovers_female_mean))))
    max(data.frame(unclass(rle(subset(complete_windowed_dataset, complete_windowed_dataset$LG=="LG1")$crossovers_female_mean))))
    
    
    
    complete_windowed_dataset$cm_mb_male_very_low<-(complete_windowed_dataset$cm_mb_male_mean<0.01)
    windowed_map$cm_mb_female_very_low<-(windowed_map$cm_mb_female<0.01)
    
    
    
    
    
    complete_windowed_dataset<-select(complete_windowed_dataset, -cm_mb_male_very_low)
    
    
    
    
    
    plot(complete_windowed_dataset$cm_mb_male, complete_windowed_dataset$pollen_biased_sum)
    
    plot(complete_windowed_dataset$cm_mb_female, complete_windowed_dataset$pollen_biased_sum)
    
    plot(complete_windowed_dataset$cm_mb_female, complete_windowed_dataset$n_BUSCOs)
    
    plot(complete_windowed_dataset$cm_mb_female, complete_windowed_dataset$sum_BUSCO_Length)
    
    View(complete_windowed_dataset %>% group_by(LG) %>% summarise(sum=sum( sum_length_RNA_TEs_all)))
    
    View(complete_windowed_dataset %>% group_by(LG) %>% summarise(mean=mean( sum_length_RNA_TEs_all)))
    
    

    
################ 4. Test predictions with linear models ################
    
    #rm(complete_windowed_dataset)
    
#    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
 #   complete_windowed_dataset<-complete_windowed_dataset%>%mutate_if(is.numeric, list(~na_if(., Inf)))
  #  complete_windowed_dataset<-complete_windowed_dataset%>%mutate_if(is.numeric, list(~na_if(., -Inf)))
   # test<-complete_windowed_dataset%>%mutate_if(is.numeric, list(~na_if(., Inf)))

################ 4.1 Inspect distributions of variables of interest ################
    
    ## Recombination rates ##
    
    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    for (i in LGs){
      print(i)
      hist(subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$cm_mb_male_mean, main=paste("cM/Mb male, LG ", i))
      hist(subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$crossovers_male_sum, main=paste("crossovers male, LG ", i))
    }

    hist(complete_windowed_dataset$cm_mb_male, breaks=100)
    hist(complete_windowed_dataset$m_crossovers, breaks=100)
    sum(as.numeric(complete_windowed_dataset$cm_mb_male==0),na.rm=T)
    
    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    for (i in LGs){
      print(i)
      hist(subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$m_crossovers, main=paste("crossovers male, LG ", i))
      
    }
    
    
    
    
    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    for (i in LGs){
      print(i)
      plot(subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$position_window,subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$cm_mb_male, main=paste("cM/Mb male, LG ", i))
        
    }
    
    
    ## Gene density ###
    
    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    for (i in LGs){
      print(i)
      hist(subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$sum_Length_Yunchen, main=paste("Gene content per window, LG ", i))
    }
    
  ################ 4.1.1 Plot correlations between variables ################
    
    
    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    for (i in LGs){
      print(i)
      plot(subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$sum_Length_Yunchen, subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$cm_mb_male, main=paste("Male recombination vs gene content, LG ", i))
    }
    
        
    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    par(mfrow=c(2,3))
    for (i in LGs){
      print(i)
      plot(subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$position_window, subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$N_RNA_TEs_all, main=paste("RNA TEs, ", i), ylim=c(0,1000), xlab="MB", ylab="N RNA TEs")
    }
    
    
    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    par(mfrow=c(2,3))
    for (i in LGs){
      print(i)
      plot(subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$position_window, subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)$sum_length_RNA_TEs_all, main=paste("RNA TEs, ", i), xlab="MB", ylab="Windowed sum RNA TEs", xlim=c(0,400))
    }
    
    
################ 4.1.2 Calculate correlations between variables ################
    
    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    for (i in LGs){
      print(i)
      LG_subset<-subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)[,c(2:64)]
      Pearson_cor<-cor(LG_subset, method = "pearson", use = "pairwise.complete.obs")
      write.csv(Pearson_cor, file=paste("8-30-21_Pearson_correlation_windowed_",i,".csv"), quote=F)
    }
    
    
    LGs<-c("LG1","LG2","LG3","LG4","LG5")  
    for (i in LGs){
      print(i)
      LG_subset<-subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)[,c(2:61)]
      Pearson_cor<-cor(LG_subset, method = "pearson", use = "pairwise.complete.obs")
      pdf(file = paste("8-10_correlations_",i,".pdf"),
          width = 4, # The width of the plot in inches
          height = 4) # The height of the plot in inches
      tmp_plot<-corrplot(Pearson_cor,method="circle", tl.pos='n')
      dev.off()
    }

    data.frame(colnames(complete_windowed_dataset))    
    
    
    
    
    
 
        for (i in LGs){
      print(i)
      LG_subset<-subset(complete_windowed_dataset, complete_windowed_dataset$LG==i)[,c(2:49,54:59,61)]
      Pearson_cor<-cor(LG_subset, method = "pearson", use = "complete.obs")
     # write.csv(Pearson_cor, file=paste("8-10-21_Pearson_correlation_windowed_",i,".csv"))
    }

    
    
    hist(complete_windowed_dataset$sum_Length_Yunchen)
    hist(complete_windowed_dataset$n_genes_Yunchen)
    plot(complete_windowed_dataset$cm_mb_male,complete_windowed_dataset$sum_Length_Yunchen)

    
    
    
    
NA_enriched<-complete_windowed_dataset[,c(50:53,60)]    
View(NA_enriched)    

leaf_expression<-subset(complete_windowed_dataset, !is.na(complete_windowed_dataset$leaf_sex_log2FoldChange_mean))
View(leaf_expression)

LGs<-c("LG1","LG2","LG3","LG4","LG5")  
for (i in LGs){
  print(i)
  LG_subset<-subset(leaf_expression, leaf_expression$LG==i)[,c(2:61)]
  Pearson_cor<-cor(LG_subset, method = "pearson", use = "pairwise.complete.obs")
 # write.csv(Pearson_cor, file=paste("8-10-21_leaf_expression_Pearson_correlation_windowed_",i,".csv"))
}



    
View(LG_subset)        
View(complete_windowed_dataset)
Pearson_cor<-cor(LG_subset, method = "pearson", use = "na.or.complete")

colnames(complete_windowed_dataset)
View(Pearson_cor)

    
    Pearson_cor<-cor(complete_windowed_dataset, method = "pearson", use = "complete.obs")
    
    rcorrAllPearson<-rcorr(as.matrix(complete_windowed_dataset,type="pearson"))
    RvalsSp<-rcorrAllSpearman$r
    PvalsSp<-rcorrAllSpearman$P
    #write.csv(PvalsSp, "All_P_spearman.csv")
    #write.csv(RvalsSp, "All_R_spearman.csv")
    
    
    
    
    
    
    
        
################ 4.2 Compare models ################
    
################ 4.2.1 Compare models for male recombination rate ################
    

    #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
    
   # View(complete_windowed_dataset)
    colnames(complete_windowed_dataset)
    
    complete_windowed_dataset<-subset(complete_windowed_dataset,!is.na(complete_windowed_dataset$cm_mb_male)) #First, remove all parts of the dataset that do not include recombination
    
      
    ################ 4.2.1.1 LG1 ################

    LG1<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG1")
    View(LG1)
    
    mean(LG1$crossovers_male_sum)
    var(LG1$crossovers_male_sum)

    hist(LG1$cm_mb_male, breaks = 30)
    
    plot(LG1$position_window , LG1$crossovers_male_sum)
    plot(LG1$n_genes_Yunchen , LG1$crossovers_male_sum)
    plot(LG1$sum_Length_Yunchen , LG1$crossovers_male_sum)
    plot(LG1$male_distortion_mean, LG1$crossovers_male_sum)
    plot(LG1$N_DNA_TEs_all, LG1$crossovers_male_sum)
    plot(LG1$N_RNA_TEs_all, LG1$crossovers_male_sum)
    plot(LG1$pollen_biased_sum, LG1$crossovers_male_sum)
    plot(LG1$expressed_leaf_sum, LG1$crossovers_male_sum)
    plot(LG1$leaf_male_bias_sum, LG1$crossovers_male_sum)
    plot(LG1$flower_male_bias_sum, LG1$crossovers_male_sum)
    plot(LG1$pollen_tube_biased_sum, LG1$crossovers_male_sum)
    
    plot(LG1$position_window , LG1$leaf_male_bias_sum)
    
    fitdistr(LG1$crossovers_male_sum, 'poisson')$loglik
    fitdistr(LG1$crossovers_male_sum, 'normal')$loglik
    fitdistr(LG1$crossovers_male_sum, 'negative binomial')$loglik
    
    LG1_distortion<-subset(LG1, !is.na(male_distortion_mean))


    male_cMMb_model1.0 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen, 
      #  family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
          family = nbinom2,
    #  ziformula=~1, 
  # control=glmmTMBControl(optimizer=optim,
                         # optArgs=list(method="BFGS")),
        #data = LG1)
        data = LG1_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(male_cMMb_model1.0)
    
      
    male_cMMb_model1.1 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen + N_DNA_TEs_all, 
        # family = poisson(link = "log"), 
       #  family = tweedie(link = "log"), 
      family = nbinom2,
    #    ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      #data = LG1)
    data = LG1_distortion)
    summary(male_cMMb_model1.1)
    #Interaction term not significant
    
    anova(male_cMMb_model1.0, male_cMMb_model1.1)
    

    male_cMMb_model1.2 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen+N_RNA_TEs_all, 
      # family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      family = nbinom2,
      #    ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
#      data = LG1)
       data = LG1_distortion)
    summary(male_cMMb_model1.2)
    anova(male_cMMb_model1.0, male_cMMb_model1.2)
    anova(male_cMMb_model1.1, male_cMMb_model1.2)
    #Interaction term not significant

    male_cMMb_model1.3 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen+N_RNA_TEs_all+male_distortion_mean, 
      # family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      family = nbinom2,
   #       ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      data = LG1_distortion)
    summary(male_cMMb_model1.3)
    anova(male_cMMb_model1.0, male_cMMb_model1.3)
    anova(male_cMMb_model1.1, male_cMMb_model1.3)
    anova(male_cMMb_model1.2, male_cMMb_model1.3)
    
    
    male_cMMb_model1.4 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen+N_RNA_TEs_all*male_distortion_mean, 
      # family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      family = nbinom2,
      #    ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      data = LG1_distortion)
    summary(male_cMMb_model1.4)
    anova(male_cMMb_model1.0, male_cMMb_model1.4)
    anova(male_cMMb_model1.1, male_cMMb_model1.4)
    anova(male_cMMb_model1.2, male_cMMb_model1.4)
    anova(male_cMMb_model1.3, male_cMMb_model1.4)

    male_cMMb_model1.5 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen+N_RNA_TEs_all*male_distortion_mean+leaf_male_bias_sum, 
      # family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      family = nbinom2,
      #    ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      data = LG1_distortion)
    summary(male_cMMb_model1.5)
    anova(male_cMMb_model1.0, male_cMMb_model1.5)
    anova(male_cMMb_model1.1, male_cMMb_model1.5)
    anova(male_cMMb_model1.2, male_cMMb_model1.5)
    anova(male_cMMb_model1.3, male_cMMb_model1.5)
    anova(male_cMMb_model1.4, male_cMMb_model1.5)
    
        



    simulationOutput <- simulateResiduals(fittedModel = male_cMMb_model1.4, plot = F, n = 1000)
    plot(simulationOutput)
    testDispersion(simulationOutput)
    testZeroInflation(simulationOutput) 
    testQuantiles(simulationOutput)

    plotResiduals(simulationOutput, LG1$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG1$n_genes_Yunchen)
    plotResiduals(simulationOutput, LG1$sum_length_DNA_TEs_all)
    plotResiduals(simulationOutput, LG1$sum_length_RNA_TEs_all)
    plotResiduals(simulationOutput, LG1$male_distortion)
    plotResiduals(simulationOutput, LG1$pollen_biased_sum)
    plotResiduals(simulationOutput, LG1_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG1$pollen_biased_sum/LG1$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG1$flower_male_bias_sum)
    plotResiduals(simulationOutput,   (LG1$flower_male_bias_sum/LG1$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG1$flower_male_bias_sum/LG1$expressed_flower_sum))
    
    
    plotResiduals(simulationOutput, LG1$expressed_flower_sum)
    
    plot(complete_windowed_dataset$male_distortion_mean, (complete_windowed_dataset$pollen_biased_sum/complete_windowed_dataset$n_genes_Yunchen))

    plot(complete_windowed_dataset$male_distortion, complete_windowed_dataset$pollen_tube_biased_sum)
    

    
    ################ 4.2.1.2 LG2 ################
    
    
    LG2<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG2")
    LG2_distortion<-subset(LG2, !is.na(male_distortion_mean))
    
    mean(LG2$crossovers_male_sum)
    var(LG2$crossovers_male_sum)
    
    hist(LG2$crossovers_male_sum, breaks = 30)
    
    plot(LG2$position_window , LG2$crossovers_male_sum)
    plot(LG2$n_genes_Yunchen , LG2$crossovers_male_sum)
    plot(LG2$sum_Length_Yunchen , LG2$crossovers_male_sum)
    plot(LG2$male_distortion_mean, LG2$crossovers_male_sum)
    plot(LG2$N_DNA_TEs_all, LG2$crossovers_male_sum)
    plot(LG2$N_RNA_TEs_all, LG2$crossovers_male_sum)
    plot(LG2$pollen_biased_sum, LG2$crossovers_male_sum)
    plot(LG2$expressed_leaf_sum, LG2$crossovers_male_sum)
    plot(LG2$leaf_male_bias_sum, LG2$crossovers_male_sum)
    plot(LG2$flower_male_bias_sum, LG2$crossovers_male_sum)
    plot(LG2$pollen_tube_biased_sum, LG2$crossovers_male_sum)
    
    
    
    
    male_cMMb_model2.0 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen, 
      # family = poisson(link = "log"), 
    #    family = tweedie(link = "log"), 
      family = nbinom2,
        ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
#      data = LG2)
      data = LG2_distortion)
    summary(male_cMMb_model2.0)
    
    male_cMMb_model2.1 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen + N_DNA_TEs_all, 
      # family = poisson(link = "log"), 
      #    family = tweedie(link = "log"), 
      family = nbinom2,
      ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
#      data = LG2)
    data = LG2_distortion)
    summary(male_cMMb_model2.1)
    anova(male_cMMb_model2.0, male_cMMb_model2.1)
    
    male_cMMb_model2.2 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen * N_DNA_TEs_all, 
      # family = poisson(link = "log"), 
      #    family = tweedie(link = "log"), 
      family = nbinom2,
      #ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
#      data = LG2)
    data = LG2_distortion)
    summary(male_cMMb_model2.2)
    anova(male_cMMb_model2.0, male_cMMb_model2.2)
    anova(male_cMMb_model2.1, male_cMMb_model2.2)

    
    male_cMMb_model2.3 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen * N_RNA_TEs_all, 
      # family = poisson(link = "log"), 
      #    family = tweedie(link = "log"), 
      family = nbinom2,
  #   ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
   #  data = LG2)
    data = LG2_distortion)

    summary(male_cMMb_model2.3)
    anova(male_cMMb_model2.0, male_cMMb_model2.3)
    anova(male_cMMb_model2.1, male_cMMb_model2.3)
    anova(male_cMMb_model2.2, male_cMMb_model2.3)
    
    
    male_cMMb_model2.4 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen + N_DNA_TEs_all + N_RNA_TEs_all + n_genes_Yunchen:N_RNA_TEs_all, 
      # family = poisson(link = "log"), 
      #    family = tweedie(link = "log"), 
      family = nbinom2,
   # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
#      data = LG2)
    data = LG2_distortion)
    summary(male_cMMb_model2.4)
    anova(male_cMMb_model2.0, male_cMMb_model2.4)
    anova(male_cMMb_model2.1, male_cMMb_model2.4)
    anova(male_cMMb_model2.2, male_cMMb_model2.4)
    anova(male_cMMb_model2.3, male_cMMb_model2.4)
    
    male_cMMb_model2.5 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen + N_RNA_TEs_all + male_distortion_mean,
      # family = poisson(link = "log"), 
      #    family = tweedie(link = "log"), 
      family = nbinom2,
    #   ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      data = LG2_distortion)
     summary(male_cMMb_model2.5)
    anova(male_cMMb_model2.0, male_cMMb_model2.5)
    anova(male_cMMb_model2.1, male_cMMb_model2.5)
    anova(male_cMMb_model2.2, male_cMMb_model2.5)
    anova(male_cMMb_model2.3, male_cMMb_model2.5)
    anova(male_cMMb_model2.4, male_cMMb_model2.5)

    
    male_cMMb_model2.6 <- glmmTMB(
      crossovers_male_sum ~ 
        N_RNA_TEs_all + male_distortion_mean,
      # family = poisson(link = "log"), 
      #    family = tweedie(link = "log"), 
      family = nbinom2,
    #     ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      data = LG2_distortion)
    


    summary(male_cMMb_model2.6)
    anova(male_cMMb_model2.0, male_cMMb_model2.6)
    anova(male_cMMb_model2.1, male_cMMb_model2.6)
    anova(male_cMMb_model2.2, male_cMMb_model2.6)
    anova(male_cMMb_model2.3, male_cMMb_model2.6)
    anova(male_cMMb_model2.4, male_cMMb_model2.6)
    anova(male_cMMb_model2.5, male_cMMb_model2.6)
    
    
    male_cMMb_model2.7 <- glmmTMB(
      crossovers_male_sum ~ 
        N_RNA_TEs_all + male_distortion_mean + pollen_tube_biased_sum,
      # family = poisson(link = "log"), 
      #    family = tweedie(link = "log"), 
      family = nbinom2,
      #     ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      data = LG2_distortion)
    


    summary(male_cMMb_model2.7)
    anova(male_cMMb_model2.0, male_cMMb_model2.7)
    anova(male_cMMb_model2.1, male_cMMb_model2.7)
    anova(male_cMMb_model2.2, male_cMMb_model2.7)
    anova(male_cMMb_model2.3, male_cMMb_model2.7)
    anova(male_cMMb_model2.4, male_cMMb_model2.7)
    anova(male_cMMb_model2.5, male_cMMb_model2.7)
    anova(male_cMMb_model2.6, male_cMMb_model2.7)
    
    crossovers_male_sum ~ 
      N_RNA_TEs_all + male_distortion_mean +  
      pollen_tube_biased_sum * flower_male_bias_sum
    
    male_cMMb_model2.8 <- glmmTMB(
      crossovers_male_sum ~ 
        N_RNA_TEs_all + male_distortion_mean + pollen_tube_biased_sum*flower_male_bias_sum,
      # family = poisson(link = "log"), 
      #    family = tweedie(link = "log"), 
      family = nbinom2,
      ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      data = LG2_distortion)
    


    summary(male_cMMb_model2.8)

    anova(male_cMMb_model2.0, male_cMMb_model2.8)
    anova(male_cMMb_model2.1, male_cMMb_model2.8)
    anova(male_cMMb_model2.2, male_cMMb_model2.8)
    anova(male_cMMb_model2.3, male_cMMb_model2.8)
    anova(male_cMMb_model2.4, male_cMMb_model2.8)
    anova(male_cMMb_model2.5, male_cMMb_model2.8)
    anova(male_cMMb_model2.6, male_cMMb_model2.8)
    anova(male_cMMb_model2.7, male_cMMb_model2.8)
    
          
    
    simulationOutput <- simulateResiduals(fittedModel = male_cMMb_model2.7, plot = F, n = 1000)
    plot(simulationOutput)
    testQuantiles(simulationOutput)




    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    plotResiduals(simulationOutput, LG2_distortion$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG2_distortion$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG2_distortion$sum_length_DNA_TEs_all)

    plotResiduals(simulationOutput, LG2_distortion$sum_length_RNA_TEs_all)

    plotResiduals(simulationOutput, LG2_distortion$male_distortion)
    plotResiduals(simulationOutput, LG2_distortion$pollen_biased_sum)
    
    plotResiduals(simulationOutput,   (LG2_distortion$pollen_biased_sum/LG2$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG2_distortion$flower_male_bias_sum)
    plotResiduals(simulationOutput, LG2_distortion$leaf_male_bias_sum)
    plotResiduals(simulationOutput, LG2_distortion$female_distortion_mean)
    plotResiduals(simulationOutput, LG2_distortion$flower_directional_bias_mean)
    
    plotResiduals(simulationOutput,   (LG2_distortion$flower_male_bias_sum/LG2_distortion$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG2_distortion$flower_male_bias_sum/LG2_distortion$expressed_flower_sum))
    
    
    

    ################ 4.2.1.3 LG3 ################
    
    
    LG3<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG3")
    LG3_distortion<-subset(LG3, !is.na(male_distortion_mean))
    
    mean(LG3$crossovers_male_sum)
    var(LG3$crossovers_male_sum)
    
    hist(LG3$crossovers_male_sum, breaks = 30)
    
    plot(LG3$position_window , LG3$crossovers_male_sum)
    plot(LG3$n_genes_Yunchen , LG3$crossovers_male_sum)
    plot(LG3$sum_Length_Yunchen , LG3$crossovers_male_sum)
    plot(LG3$male_distortion_mean, LG3$crossovers_male_sum)
    plot(LG3$N_DNA_TEs_all, LG3$crossovers_male_sum)
    plot(LG3$N_RNA_TEs_all, LG3$crossovers_male_sum)
    plot(LG3$pollen_biased_sum, LG3$crossovers_male_sum)
    plot(LG3$expressed_leaf_sum, LG3$crossovers_male_sum)
    plot(LG3$leaf_male_bias_sum, LG3$crossovers_male_sum)
    plot(LG3$flower_male_bias_sum, LG3$crossovers_male_sum)
    plot(LG3$pollen_biased_sum, LG3$crossovers_male_sum)
    plot(LG3$pollen_tube_biased_sum, LG3$crossovers_male_sum)
    
    
    
    
    male_cMMb_model3.0 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen, 
      # family = poisson(link = "log"), 
     #     family = tweedie(link = "log"), 
      family = nbinom2,
    #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
        #    data = LG3)
      data = LG3_distortion)

    summary(male_cMMb_model3.0)
    
    
    
    male_cMMb_model3.1 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen + N_RNA_TEs_all, 
      # family = poisson(link = "log"), 
      #     family = tweedie(link = "log"), 
      family = nbinom2,
       # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    #  data = LG3)
    data = LG3_distortion)
    
    summary(male_cMMb_model3.1)
    
    anova(male_cMMb_model3.0,male_cMMb_model3.1)
    
    
    
    male_cMMb_model3.2 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen + male_distortion_mean, 
      # family = poisson(link = "log"), 
      #     family = tweedie(link = "log"), 
      family = nbinom2,
     #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      #data = LG3)
    data = LG3_distortion)
    
    summary(male_cMMb_model3.2)
    
    anova(male_cMMb_model3.0,male_cMMb_model3.2)
    anova(male_cMMb_model3.1,male_cMMb_model3.2)
    
    
    male_cMMb_model3.3 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen * male_distortion_mean, 
      # family = poisson(link = "log"), 
      #     family = tweedie(link = "log"), 
      family = nbinom2,
      #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      #data = LG3)
      data = LG3_distortion)
    
    summary(male_cMMb_model3.3)
    
    anova(male_cMMb_model3.0,male_cMMb_model3.3)
    anova(male_cMMb_model3.1,male_cMMb_model3.3)
    anova(male_cMMb_model3.2,male_cMMb_model3.3)
    
    male_cMMb_model3.4 <- glmmTMB(
      crossovers_male_sum ~ 
        n_genes_Yunchen * male_distortion_mean + N_RNA_TEs_all, 
      # family = poisson(link = "log"), 
      #     family = tweedie(link = "log"), 
      family = nbinom2,
      #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      #data = LG3)
      data = LG3_distortion)
    
    summary(male_cMMb_model3.4)
    
    female_distortion_mean
    
    anova(male_cMMb_model3.0,male_cMMb_model3.4)
    anova(male_cMMb_model3.1,male_cMMb_model3.4)
    anova(male_cMMb_model3.2,male_cMMb_model3.4)
    anova(male_cMMb_model3.3,male_cMMb_model3.4)
    
    
 
    simulationOutput <- simulateResiduals(fittedModel = male_cMMb_model3.4, plot = F, n = 1000)
    plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    plotResiduals(simulationOutput, LG3$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG3$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG3$sum_length_DNA_TEs_all)

    plotResiduals(simulationOutput, LG3$sum_length_RNA_TEs_all)

    plotResiduals(simulationOutput, LG3$expressed_leaf_sum)  
    plotResiduals(simulationOutput, LG3$expressed_flower_sum)  
    plotResiduals(simulationOutput, LG3$male_distortion_mean)  
    plotResiduals(simulationOutput, LG3$leaf_male_bias_mean)  
    plotResiduals(simulationOutput, LG3$leaf_male_bias_sum)  
    plotResiduals(simulationOutput, LG3$flower_male_bias_sum)
    plotResiduals(simulationOutput, LG3$pollen_tube_biased_sum)  
    plotResiduals(simulationOutput, LG3_distortion$position_window)  
    
    View(LG3)
    

    
    ################ 4.2.1.4 LG4 ################
    
    
    LG4<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG4")
    LG4_distortion<-subset(LG4, !is.na(male_distortion_mean))
    
    mean(LG4$crossovers_male_sum)
    var(LG4$crossovers_male_sum)
    
    hist(LG4$crossovers_male_sum, breaks = 30)
    
    plot(LG4$position_window , LG4$crossovers_male_sum)
    plot(LG4$n_genes_Yunchen , LG4$crossovers_male_sum)
    plot(LG4$sum_Length_Yunchen , LG4$crossovers_male_sum)
    plot(LG4$male_distortion_mean, LG4$crossovers_male_sum)
    plot(LG4$N_DNA_TEs_all, LG4$crossovers_male_sum)
    plot(LG4$N_RNA_TEs_all, LG4$crossovers_male_sum)
    plot(LG4$pollen_biased_sum, LG4$crossovers_male_sum)
    plot(LG4$expressed_leaf_sum, LG4$crossovers_male_sum)
    plot(LG4$leaf_male_bias_sum, LG4$crossovers_male_sum)
    plot(LG4$flower_male_bias_sum, LG4$crossovers_male_sum)
    plot(LG4$pollen_biased_sum, LG4$crossovers_male_sum)
    plot(LG4$pollen_tube_biased_sum, LG4$crossovers_male_sum)
    plot(LG4$pollen_tube_biased_sum, LG4$crossovers_female_sum)
    
    
    
    male_cMMb_model4.0 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window, 
      # family = poisson(link = "log"), 
           family = tweedie(link = "log"), 
      #family = nbinom2,
      #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
     #     data = LG4)
      data = LG4_distortion)
    summary(male_cMMb_model4.0)
    #residuals plots suggest m distortion and leaf male bias important
 
    male_cMMb_model4.1 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window * n_genes_Yunchen, 
      # family = poisson(link = "log"), 
      family = tweedie(link = "log"), 
      #family = nbinom2,
       # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
       #   data = LG4)
      data = LG4_distortion)
    summary(male_cMMb_model4.1)
    #residuals plots suggest leaf male bias and pollen tube bias important; male distortion in reduced dataset
    anova(male_cMMb_model4.0,male_cMMb_model4.1)
    
    
    male_cMMb_model4.2 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window * n_genes_Yunchen + pollen_biased_sum, 
      # family = poisson(link = "log"), 
      family = tweedie(link = "log"), 
      #family = nbinom2,
     #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
       #data = LG4) #ADD ZERO INFLATION FOR FULL LG
      data = LG4_distortion) #REMOVE ZERO INFLATION FOR DISTORTION SUBSET
    summary(male_cMMb_model4.2)
    #residuals plots suggest leaf male bias and pollen tube bias important; male distortion in reduced dataset
    anova(male_cMMb_model4.0,male_cMMb_model4.2)
    anova(male_cMMb_model4.1,male_cMMb_model4.2)
    
    male_cMMb_model4.3 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window * n_genes_Yunchen + n_genes_Yunchen:male_distortion_mean + pollen_biased_sum, 
      # family = poisson(link = "log"), 
      family = tweedie(link = "log"), 
      #family = nbinom2,
      #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      #data = LG4)
    data = LG4_distortion)
    summary(male_cMMb_model4.3)
    #residuals plots suggest leaf male bias, pollen tube bias important
    anova(male_cMMb_model4.0,male_cMMb_model4.3)
    anova(male_cMMb_model4.1,male_cMMb_model4.3)
    anova(male_cMMb_model4.2,male_cMMb_model4.3)
    
    
    male_cMMb_model4.4 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window * n_genes_Yunchen + n_genes_Yunchen:male_distortion_mean + pollen_biased_sum *pollen_tube_biased_sum + leaf_male_bias_sum, 
      # family = poisson(link = "log"), 
      family = tweedie(link = "log"), 
      #family = nbinom2,
      #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      #data = LG4)
      data = LG4_distortion)
    summary(male_cMMb_model4.4)
    #residuals plots suggest leaf male bias, pollen tube bias important
    anova(male_cMMb_model4.0,male_cMMb_model4.4)
    anova(male_cMMb_model4.1,male_cMMb_model4.4)
    anova(male_cMMb_model4.2,male_cMMb_model4.4)
    anova(male_cMMb_model4.3,male_cMMb_model4.4)
    
    
    
    
    simulationOutput <- simulateResiduals(fittedModel = male_cMMb_model4.4, plot = F, n = 1000)
    plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    plotResiduals(simulationOutput, LG4$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG4$n_genes_Yunchen)
    plotResiduals(simulationOutput, LG4_distortion$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG4$sum_length_DNA_TEs_all)

    plotResiduals(simulationOutput, LG4$sum_length_RNA_TEs_all)
    plotResiduals(simulationOutput, LG4_distortion$sum_length_RNA_TEs_all)

    plotResiduals(simulationOutput, LG4$expressed_leaf_sum)  
    plotResiduals(simulationOutput, LG4$expressed_flower_sum)  
    plotResiduals(simulationOutput, LG4$male_distortion_mean)  
    plotResiduals(simulationOutput, LG4$leaf_male_bias_mean)  
    plotResiduals(simulationOutput, LG4$leaf_male_bias_sum)  
    plotResiduals(simulationOutput, LG4$flower_male_bias_sum)
    plotResiduals(simulationOutput, LG4$pollen_tube_biased_sum)  
    plotResiduals(simulationOutput, LG4$position_window)  
    
    plotResiduals(simulationOutput, LG4_distortion$expressed_leaf_sum)  
    plotResiduals(simulationOutput, LG4_distortion$expressed_flower_sum)  
    plotResiduals(simulationOutput, LG4_distortion$male_distortion_mean)  
    plotResiduals(simulationOutput, LG4_distortion$leaf_male_bias_mean)  
    plotResiduals(simulationOutput, LG4_distortion$leaf_male_bias_sum)  
    plotResiduals(simulationOutput, LG4_distortion$flower_male_bias_sum)
    plotResiduals(simulationOutput, LG4_distortion$pollen_tube_biased_sum)  
    plotResiduals(simulationOutput, LG4_distortion$position_window)  
    
    
    
    
    ################ 4.2.1.5 LG5 ################
    
    
    LG5<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG5")
    LG5_distortion<-subset(LG5, !is.na(male_distortion_mean))
    
    mean(LG5$crossovers_male_sum)
    var(LG5$crossovers_male_sum)
    
    hist(LG5$crossovers_male_sum, breaks = 30)
    
    plot(LG5$position_window , LG5$crossovers_male_sum)
    plot(LG5$n_genes_Yunchen , LG5$crossovers_male_sum)
    plot(LG5$sum_Length_Yunchen , LG5$crossovers_male_sum)
    plot(LG5$male_distortion_mean, LG5$crossovers_male_sum)
    plot(LG5$N_DNA_TEs_all, LG5$crossovers_male_sum)
    plot(LG5$N_RNA_TEs_all, LG5$crossovers_male_sum)
    plot(LG5$pollen_biased_sum, LG5$crossovers_male_sum)
    plot(LG5$expressed_leaf_sum, LG5$crossovers_male_sum)
    plot(LG5$leaf_male_bias_sum, LG5$crossovers_male_sum)
    plot(LG5$flower_male_bias_sum, LG5$crossovers_male_sum)
    plot(LG5$pollen_biased_sum, LG5$crossovers_male_sum)
    plot(LG5$pollen_tube_biased_sum, LG5$crossovers_male_sum)
    plot(LG5$pollen_tube_biased_sum, LG5$crossovers_female_sum)
    
    
    
    male_cMMb_model5.0 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window*n_genes_Yunchen, 
    #   family = poisson(link = "log"), 
     # family = tweedie(link = "log"), 
      family = nbinom2,
       # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    #       data = LG5)
     data = LG5_distortion)
    summary(male_cMMb_model5.0)
    #residuals plots suggest m distortion and leaf male bias important
    
    male_cMMb_model5.1 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window*n_genes_Yunchen + pollen_tube_biased_sum, 
      #   family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    #  data = LG5)
    data = LG5_distortion)
    summary(male_cMMb_model5.1)
    #residuals plots suggest m distortion and leaf male bias important
    anova(male_cMMb_model5.0,male_cMMb_model5.1)
    
    
    male_cMMb_model5.2 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window*n_genes_Yunchen + leaf_male_bias_sum + pollen_tube_biased_sum, 
      #   family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
       #     data = LG5)
      data = LG5_distortion)
    summary(male_cMMb_model5.2)
    #residuals plots suggest m distortion and leaf male bias important
    anova(male_cMMb_model5.0,male_cMMb_model5.2)
    anova(male_cMMb_model5.1,male_cMMb_model5.2)
    
    
    male_cMMb_model5.3 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window*n_genes_Yunchen*male_distortion_mean + pollen_tube_biased_sum, 
      #   family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      #     data = LG5)
      data = LG5_distortion)
    summary(male_cMMb_model5.3)
    anova(male_cMMb_model5.0,male_cMMb_model5.3)
    anova(male_cMMb_model5.1,male_cMMb_model5.3)
    anova(male_cMMb_model5.2,male_cMMb_model5.3)
    
    
    
    male_cMMb_model5.4 <- glmmTMB(
      crossovers_male_sum ~ 
        position_window*male_distortion_mean + pollen_tube_biased_sum, 
      #   family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      #     data = LG5)
      data = LG5_distortion)
    summary(male_cMMb_model5.4)
    anova(male_cMMb_model5.0,male_cMMb_model5.4)
    anova(male_cMMb_model5.1,male_cMMb_model5.4)
    anova(male_cMMb_model5.2,male_cMMb_model5.4)
    anova(male_cMMb_model5.3,male_cMMb_model5.4)
    
    
    
    
    simulationOutput <- simulateResiduals(fittedModel = male_cMMb_model5.4, plot = F, n = 1000)
    plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    plotResiduals(simulationOutput, LG5$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG5$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG5$n_DNA_TEs_all)

    plotResiduals(simulationOutput, LG5$n_RNA_TEs_all)

    plotResiduals(simulationOutput, LG5$expressed_leaf_sum)  
    plotResiduals(simulationOutput, LG5$expressed_flower_sum)  
    plotResiduals(simulationOutput, LG5$male_distortion_mean)  
    plotResiduals(simulationOutput, LG5$leaf_male_bias_mean)  
    plotResiduals(simulationOutput, LG5$leaf_male_bias_sum)  
    plotResiduals(simulationOutput, LG5$flower_male_bias_sum)
    plotResiduals(simulationOutput, LG5$pollen_tube_biased_sum)  
    plotResiduals(simulationOutput, LG5$pollen_biased_sum)  
    plotResiduals(simulationOutput, LG5$position_window)  
    
    plotResiduals(simulationOutput, LG5_distortion$n_genes_Yunchen)
    plotResiduals(simulationOutput, LG5_distortion$sum_length_DNA_TEs_all)
    plotResiduals(simulationOutput, LG5_distortion$sum_length_RNA_TEs_all)
    plotResiduals(simulationOutput, LG5_distortion$expressed_leaf_sum)  
    plotResiduals(simulationOutput, LG5_distortion$expressed_flower_sum)  
    plotResiduals(simulationOutput, LG5_distortion$male_distortion_mean)  
    plotResiduals(simulationOutput, LG5_distortion$leaf_male_bias_mean)  
    plotResiduals(simulationOutput, LG5_distortion$leaf_male_bias_sum)  
    plotResiduals(simulationOutput, LG5_distortion$flower_male_bias_sum)
    plotResiduals(simulationOutput, LG5_distortion$pollen_tube_biased_sum)  
    plotResiduals(simulationOutput, LG5_distortion$position_window)  
    plotResiduals(simulationOutput, LG5_distortion$crossovers_female_sum)  
    
    
    
    
    ################ 4.2.2 Compare models for female recombination rate ################
    
    #Mostly, the distribution of cM/Mb is vehemently not normal, so I will need to find the appropriate distribution - probably Tweedie distributions in glmmTMB, since those are suitable for "skewed continuous responses with a spike at zero"
    
    #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
    
    # View(complete_windowed_dataset)
    colnames(complete_windowed_dataset)
    

    #Models including all LGs with an effect for LG don't converge, which maybe makes sense because 1) they get really overparameterized really fast and 2) different things are probably happening on the different chromosomes
    
    ################ 4.2.2.1 LG1 ################
    
    LG1<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG1")
    View(LG1)
    
    mean(LG1$crossovers_female_sum)
    var(LG1$crossovers_female_sum)
    
    hist(LG1$crossovers_female_sum, breaks = 30)
    
    plot(LG1$position_window , LG1$crossovers_female_sum)
    plot(LG1$n_genes_Yunchen , LG1$crossovers_female_sum)
    plot(LG1$sum_Length_Yunchen , LG1$crossovers_female_sum)
    plot(LG1$female_distortion_mean, LG1$crossovers_female_sum)
    plot(LG1$N_DNA_TEs_all, LG1$crossovers_female_sum)
    plot(LG1$N_RNA_TEs_all, LG1$crossovers_female_sum)
    plot(LG1$pollen_biased_sum, LG1$crossovers_female_sum)
    plot(LG1$expressed_leaf_sum, LG1$crossovers_female_sum)
    plot(LG1$leaf_female_bias_sum, LG1$crossovers_female_sum)
    plot(LG1$flower_female_bias_sum, LG1$crossovers_female_sum)
    plot(LG1$pollen_tube_biased_sum, LG1$crossovers_female_sum)
    plot(LG1$position_window , LG1$pollen_tube_biased_sum)
    
    
    fitdistr(LG1$crossovers_female_sum, 'poisson')$loglik
    fitdistr(LG1$crossovers_female_sum, 'normal')$loglik
    fitdistr(LG1$crossovers_female_sum, 'negative binomial')$loglik
    
    LG1_distortion<-subset(LG1, !is.na(female_distortion_mean))
    
    
    female_cMMb_model1.0 <- glmmTMB(
      crossovers_female_sum ~ 
        n_genes_Yunchen, 
      #  family = poisson(link = "log"), 
       # family = tweedie(link = "log"), 
      family = nbinom2,
     #   ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
     # data = LG1)
     data = LG1_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model1.0)
    
 
    female_cMMb_model1.1 <- glmmTMB(
      crossovers_female_sum ~ 
        n_genes_Yunchen * N_RNA_TEs_all, 
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
#       ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
   #  data = LG1)
     data = LG1_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model1.1)
    anova(female_cMMb_model1.0,female_cMMb_model1.1)

    
    female_cMMb_model1.2 <- glmmTMB(
      crossovers_female_sum ~ 
        n_genes_Yunchen + N_RNA_TEs_all * pollen_tube_biased_sum, 
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      #       ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
     #   data = LG1)
      data = LG1_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model1.2)
    anova(female_cMMb_model1.0,female_cMMb_model1.2)
    anova(female_cMMb_model1.1,female_cMMb_model1.2)
    
    
    
    
    female_cMMb_model1.3 <- glmmTMB(
      crossovers_female_sum ~ 
        n_genes_Yunchen+female_distortion_mean + N_RNA_TEs_all * pollen_tube_biased_sum , 
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      #       ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
       #  data = LG1)
      data = LG1_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model1.3)
    anova(female_cMMb_model1.0,female_cMMb_model1.3)
    anova(female_cMMb_model1.1,female_cMMb_model1.3)
    anova(female_cMMb_model1.2,female_cMMb_model1.3)
    
    
    
    
    
    
    
        
    simulationOutput <- simulateResiduals(fittedModel = female_cMMb_model1.3, plot = F, n = 1000)
    plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    plotResiduals(simulationOutput, LG1$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG1$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG1$sum_length_DNA_TEs_all)

    plotResiduals(simulationOutput, LG1$sum_length_RNA_TEs_all)

    plotResiduals(simulationOutput, LG1$female_distortion)
    plotResiduals(simulationOutput, LG1$pollen_biased_sum)
    plotResiduals(simulationOutput, LG1_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG1$pollen_biased_sum/LG1$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG1$flower_female_bias_sum)
    plotResiduals(simulationOutput, LG1$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG1$expressed_flower_sum)
    
    
    plotResiduals(simulationOutput,   (LG1$flower_female_bias_sum/LG1$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG1$flower_female_bias_sum/LG1$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG1_distortion$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG1_distortion$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG1_distortion$n_DNA_TEs_all)

    plotResiduals(simulationOutput, LG1_distortion$n_RNA_TEs_all)

    plotResiduals(simulationOutput, LG1_distortion$female_distortion)
    plotResiduals(simulationOutput, LG1_distortion$pollen_biased_sum)
    plotResiduals(simulationOutput, LG1_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG1_distortion$pollen_biased_sum/LG1_distortion$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG1_distortion$flower_female_bias_sum)
    plotResiduals(simulationOutput,   (LG1_distortion$flower_female_bias_sum/LG1_distortion$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG1_distortion$flower_female_bias_sum/LG1_distortion$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG1_distortion$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG1_distortion$expressed_flower_sum)
    
    plot(complete_windowed_dataset$male_distortion_mean, (complete_windowed_dataset$pollen_biased_sum/complete_windowed_dataset$n_genes_Yunchen))
    
    plot(complete_windowed_dataset$male_distortion, complete_windowed_dataset$pollen_tube_biased_sum)
    
    ################ 4.2.2.2 LG2 ################
    
    LG2<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG2")
    View(LG2)
    
    mean(LG2$crossovers_female_sum)
    var(LG2$crossovers_female_sum)
    
    hist(LG2$crossovers_female_sum, breaks = 30)
    hist(LG2$female_distortion_mean, breaks = 30)
    
    plot(LG2$position_window , LG2$crossovers_female_sum)
    plot(LG2$n_genes_Yunchen , LG2$crossovers_female_sum)
    plot(LG2$sum_Length_Yunchen , LG2$crossovers_female_sum)
    plot(LG2$female_distortion_mean, LG2$crossovers_female_sum)
    plot(LG2$N_DNA_TEs_all, LG2$crossovers_female_sum)
    plot(LG2$N_RNA_TEs_all, LG2$crossovers_female_sum)
    plot(LG2$pollen_biased_sum, LG2$crossovers_female_sum)
    plot(LG2$expressed_leaf_sum, LG2$crossovers_female_sum)
    plot(LG2$leaf_female_bias_sum, LG2$crossovers_female_sum)
    plot(LG2$flower_female_bias_sum, LG2$crossovers_female_sum)
    plot(LG2$pollen_tube_biased_sum, LG2$crossovers_female_sum)
    plot(LG2$position_window , LG2$pollen_tube_biased_sum)
    plot(LG2$position_window , LG2$female_distortion_mean)
    
    
    fitdistr(LG2$crossovers_female_sum, 'poisson')$loglik
    fitdistr(LG2$crossovers_female_sum, 'normal')$loglik
    fitdistr(LG2$crossovers_female_sum, 'negative binomial')$loglik
    
    LG2_distortion<-subset(LG2, !is.na(female_distortion_mean))
    
    
    female_cMMb_model2.0 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window * n_genes_Yunchen,
    #    family = poisson(link = "log"), 
       #family = tweedie(link = "log"), 
      family = nbinom2,
         ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
       data = LG2)
    #  data = LG2_distortion)
    #negative binomial slightly better with zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model2.0)
    
    
    female_cMMb_model2.1 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window * n_genes_Yunchen,
      #    family = poisson(link = "log"), 
      #family = tweedie(link = "log"), 
      family = nbinom2,
    #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    #  data = LG2)
    data = LG2_distortion)
   #run on both full and distortion-only dataset
    summary(female_cMMb_model2.1)
    anova(female_cMMb_model2.0,female_cMMb_model2.1)
    
    colnames(LG2)
    
    
    
    
    female_cMMb_model2.2 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window * n_genes_Yunchen + female_distortion_mean,
      #    family = poisson(link = "log"), 
      #family = tweedie(link = "log"), 
      family = nbinom2,
      #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    #  data = LG2)
    data = LG2_distortion)
    #run on both full and distortion-only dataset
    summary(female_cMMb_model2.2)
    anova(female_cMMb_model2.1,female_cMMb_model2.2)
    
    
    
    simulationOutput <- simulateResiduals(fittedModel = female_cMMb_model2.2, plot = F, n = 1000)
    plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    plotResiduals(simulationOutput, LG2$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG2$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG2$sum_length_DNA_TEs_all)

    plotResiduals(simulationOutput, LG2$sum_length_RNA_TEs_all)

    plotResiduals(simulationOutput, LG2$female_distortion)
    plotResiduals(simulationOutput, LG2$pollen_biased_sum)
    plotResiduals(simulationOutput, LG2_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG2$pollen_biased_sum/LG2$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG2$flower_female_bias_sum)
    plotResiduals(simulationOutput, LG2$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG2$expressed_flower_sum)
    
    
    plotResiduals(simulationOutput,   (LG2$flower_female_bias_sum/LG2$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG2$flower_female_bias_sum/LG2$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG2_distortion$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG2_distortion$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG2_distortion$n_DNA_TEs_all)

    plotResiduals(simulationOutput, LG2_distortion$n_RNA_TEs_all)

    plotResiduals(simulationOutput, LG2_distortion$female_distortion)
    plotResiduals(simulationOutput, LG2_distortion$pollen_biased_sum)
    plotResiduals(simulationOutput, LG2_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG2_distortion$pollen_biased_sum/LG2_distortion$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG2_distortion$flower_female_bias_sum)
    plotResiduals(simulationOutput,   (LG2_distortion$flower_female_bias_sum/LG2_distortion$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG2_distortion$flower_female_bias_sum/LG2_distortion$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG2_distortion$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG2_distortion$expressed_flower_sum)

    
    ################ 4.2.2.3 LG3 ################
    
    LG3<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG3")
    #View(LG3)
    
    mean(LG3$crossovers_female_sum)
    var(LG3$crossovers_female_sum)
    
    hist(LG3$crossovers_female_sum, breaks = 30)
    hist(LG3$female_distortion_mean, breaks = 30)
    
    plot(LG3$position_window , LG3$crossovers_female_sum)
    plot(LG3$n_genes_Yunchen , LG3$crossovers_female_sum)
    plot(LG3$sum_Length_Yunchen , LG3$crossovers_female_sum)
    plot(LG3$female_distortion_mean, LG3$crossovers_female_sum)
    plot(LG3$N_DNA_TEs_all, LG3$crossovers_female_sum)
    plot(LG3$N_RNA_TEs_all, LG3$crossovers_female_sum)
    plot(LG3$pollen_biased_sum, LG3$crossovers_female_sum)
    plot(LG3$expressed_leaf_sum, LG3$crossovers_female_sum)
    plot(LG3$leaf_female_bias_sum, LG3$crossovers_female_sum)
    plot(LG3$flower_female_bias_sum, LG3$crossovers_female_sum)
    plot(LG3$pollen_tube_biased_sum, LG3$crossovers_female_sum)
    plot(LG3$position_window , LG3$pollen_tube_biased_sum)
    
    
    fitdistr(LG3$crossovers_female_sum, 'poisson')$loglik
    fitdistr(LG3$crossovers_female_sum, 'normal')$loglik
    fitdistr(LG3$crossovers_female_sum, 'negative binomial')$loglik
    
    LG3_distortion<-subset(LG3, !is.na(female_distortion_mean))
    
        
    female_cMMb_model3.0 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window*n_genes_Yunchen,
    #  n_genes_Yunchen,
       #  family = poisson(link = "log"), 
     # family = tweedie(link = "log"), 
      family = nbinom2,
    #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
     # data = LG3)
      data = LG3_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model3.0)
    
    
    
    female_cMMb_model3.1 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window*n_genes_Yunchen + female_distortion_mean,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      #  ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    #   data = LG3)
      data = LG3_distortion)
    #negative binomial slightly better without zero-inflation
    summary(female_cMMb_model3.1)
    anova(female_cMMb_model3.0,female_cMMb_model3.1)
    
    
    
    
    
    simulationOutput <- simulateResiduals(fittedModel = female_cMMb_model3.1, plot = F, n = 1000)
    plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    plotResiduals(simulationOutput, LG3$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG3$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG3$sum_length_DNA_TEs_all)

    plotResiduals(simulationOutput, LG3$sum_length_RNA_TEs_all)

    plotResiduals(simulationOutput, LG3$female_distortion)
    plotResiduals(simulationOutput, LG3$pollen_biased_sum)
    plotResiduals(simulationOutput, LG3_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG3$pollen_biased_sum/LG3$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG3$flower_female_bias_sum)
    plotResiduals(simulationOutput, LG3$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG3$expressed_flower_sum)
    
    
    plotResiduals(simulationOutput,   (LG3$flower_female_bias_sum/LG3$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG3$flower_female_bias_sum/LG3$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG3_distortion$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG3_distortion$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG3_distortion$n_DNA_TEs_all)

    plotResiduals(simulationOutput, LG3_distortion$n_RNA_TEs_all)

    plotResiduals(simulationOutput, LG3_distortion$female_distortion)
    plotResiduals(simulationOutput, LG3_distortion$pollen_biased_sum)
    plotResiduals(simulationOutput, LG3_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG3_distortion$pollen_biased_sum/LG3_distortion$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG3_distortion$flower_female_bias_sum)
    plotResiduals(simulationOutput,   (LG3_distortion$flower_female_bias_sum/LG3_distortion$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG3_distortion$flower_female_bias_sum/LG3_distortion$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG3_distortion$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG3_distortion$expressed_flower_sum)

    
    ################ 4.2.2.4 LG4 ################
    
    LG4<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG4")
    #View(LG4)
    
    mean(LG4$crossovers_female_sum)
    var(LG4$crossovers_female_sum)
    
    hist(LG4$crossovers_female_sum, breaks = 30)
    hist(LG4$female_distortion_mean, breaks = 30)
    
    plot(LG4$n_BUSCOs , LG4$crossovers_female_sum)
    plot(LG4$position_window , LG4$crossovers_female_sum)
    plot(LG4$n_genes_Yunchen , LG4$crossovers_female_sum)
    plot(LG4$sum_Length_Yunchen , LG4$crossovers_female_sum)
    plot(LG4$female_distortion_mean, LG4$crossovers_female_sum)
    plot(LG4$N_DNA_TEs_all, LG4$crossovers_female_sum)
    plot(LG4$N_RNA_TEs_all, LG4$crossovers_female_sum)
    plot(LG4$pollen_biased_sum, LG4$crossovers_female_sum)
    plot(LG4$expressed_leaf_sum, LG4$crossovers_female_sum)
    plot(LG4$leaf_female_bias_sum, LG4$crossovers_female_sum)
    plot(LG4$flower_female_bias_sum, LG4$crossovers_female_sum)
    plot(LG4$pollen_tube_biased_sum, LG4$crossovers_female_sum)
    plot(LG4$position_window , LG4$pollen_tube_biased_sum)
    
    
    fitdistr(LG4$crossovers_female_sum, 'poisson')$loglik
    fitdistr(LG4$crossovers_female_sum, 'normal')$loglik
    fitdistr(LG4$crossovers_female_sum, 'negative binomial')$loglik
    
    LG4_distortion<-subset(LG4, !is.na(female_distortion_mean))
    
  
    female_cMMb_model4.0 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window+n_genes_Yunchen,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
     #  family = tweedie(link = "log"), 
      family = nbinom2,
       # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    #   data = LG4)
      data = LG4_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model4.0)
    
    
    female_cMMb_model4.1 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window*flower_female_bias_sum+n_genes_Yunchen,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      #   data = LG4)
      data = LG4_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model4.1)
    anova(female_cMMb_model4.0, female_cMMb_model4.1)
    
    
    female_cMMb_model4.2 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window*flower_female_bias_sum+pollen_tube_biased_sum+n_genes_Yunchen,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
   #    data = LG4)
      data = LG4_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model4.2)
    anova(female_cMMb_model4.0, female_cMMb_model4.2)
    anova(female_cMMb_model4.1, female_cMMb_model4.2)
    
    female_cMMb_model4.3 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window*flower_female_bias_sum*pollen_tube_biased_sum+n_genes_Yunchen,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    #  data = LG4)
      data = LG4_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model4.3)
    anova(female_cMMb_model4.0, female_cMMb_model4.3)
    anova(female_cMMb_model4.1, female_cMMb_model4.3)
    anova(female_cMMb_model4.2, female_cMMb_model4.3)
    
    
    female_cMMb_model4.4 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window*flower_female_bias_sum+pollen_tube_biased_sum+n_genes_Yunchen + female_distortion_mean,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    #    data = LG4)
      data = LG4_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model4.4)
    anova(female_cMMb_model4.0, female_cMMb_model4.4)
    anova(female_cMMb_model4.1, female_cMMb_model4.4)
    anova(female_cMMb_model4.2, female_cMMb_model4.4)
    anova(female_cMMb_model4.3, female_cMMb_model4.4)
    
    
    
    
    
    simulationOutput <- simulateResiduals(fittedModel = female_cMMb_model4.4, plot = F, n = 1000)
    plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    plotResiduals(simulationOutput, LG4$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG4$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG4$sum_length_DNA_TEs_all)

    plotResiduals(simulationOutput, LG4$sum_length_RNA_TEs_all)

    plotResiduals(simulationOutput, LG4$female_distortion)
    plotResiduals(simulationOutput, LG4$pollen_biased_sum)
    plotResiduals(simulationOutput, LG4_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG4$pollen_biased_sum/LG4$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG4$flower_female_bias_sum)
    plotResiduals(simulationOutput, LG4$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG4$expressed_flower_sum)
    
    
    plotResiduals(simulationOutput,   (LG4$flower_female_bias_sum/LG4$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG4$flower_female_bias_sum/LG4$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG4_distortion$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG4_distortion$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG4_distortion$n_DNA_TEs_all)

    plotResiduals(simulationOutput, LG4_distortion$n_RNA_TEs_all)

    plotResiduals(simulationOutput, LG4_distortion$female_distortion)
    plotResiduals(simulationOutput, LG4_distortion$pollen_biased_sum)
    plotResiduals(simulationOutput, LG4_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG4_distortion$pollen_biased_sum/LG4_distortion$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG4_distortion$flower_female_bias_sum)
    plotResiduals(simulationOutput,   (LG4_distortion$flower_female_bias_sum/LG4_distortion$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG4_distortion$flower_female_bias_sum/LG4_distortion$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG4_distortion$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG4_distortion$expressed_flower_sum)
    
      
    ################ 4.2.2.5 LG5 ################
    
    LG5<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG5")
    #View(LG5)
    
    mean(LG5$crossovers_female_sum)
    var(LG5$crossovers_female_sum)
    
    hist(LG5$crossovers_female_sum, breaks = 30)
    hist(LG5$female_distortion_mean, breaks = 30)
    
    plot(LG5$position_window , LG5$crossovers_female_sum)
    plot(LG5$n_genes_Yunchen , LG5$crossovers_female_sum)
    plot(LG5$sum_Length_Yunchen , LG5$crossovers_female_sum)
    plot(LG5$female_distortion_mean, LG5$crossovers_female_sum)
    plot(LG5$N_DNA_TEs_all, LG5$crossovers_female_sum)
    plot(LG5$N_RNA_TEs_all, LG5$crossovers_female_sum)
    plot(LG5$pollen_biased_sum, LG5$crossovers_female_sum)
    plot(LG5$expressed_leaf_sum, LG5$crossovers_female_sum)
    plot(LG5$leaf_female_bias_sum, LG5$crossovers_female_sum)
    plot(LG5$flower_female_bias_sum, LG5$crossovers_female_sum)
    plot(LG5$pollen_tube_biased_sum, LG5$crossovers_female_sum)
    plot(LG5$position_window , LG5$pollen_tube_biased_sum)
    
    
    fitdistr(LG5$crossovers_female_sum, 'poisson')$loglik
    fitdistr(LG5$crossovers_female_sum, 'normal')$loglik
    fitdistr(LG5$crossovers_female_sum, 'negative binomial')$loglik
    
    LG5_distortion<-subset(LG5, !is.na(female_distortion_mean))
    
    
    female_cMMb_model5.0 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window * n_genes_Yunchen,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
       # family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
         data = LG5)
     # data = LG5_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model5.0)
    
    
    
    female_cMMb_model5.1 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window * n_genes_Yunchen+N_RNA_TEs_all,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
    
       data = LG5)
   #  data = LG5_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model5.1)
    anova(female_cMMb_model5.0, female_cMMb_model5.1)
    
    
    female_cMMb_model5.2 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window * n_genes_Yunchen+position_window*expressed_leaf_sum,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
      data = LG5)
  #    data = LG5_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model5.2)
    anova(female_cMMb_model5.0, female_cMMb_model5.2)
    anova(female_cMMb_model5.1, female_cMMb_model5.2)
    
    
    female_cMMb_model5.3 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window * n_genes_Yunchen+position_window*expressed_leaf_sum +expressed_flower_sum,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
       data = LG5)
    #  data = LG5_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model5.3)
    anova(female_cMMb_model5.0, female_cMMb_model5.3)
    anova(female_cMMb_model5.1, female_cMMb_model5.3)
    anova(female_cMMb_model5.2, female_cMMb_model5.3)
    
    
    female_cMMb_model5.4 <- glmmTMB(
      crossovers_female_sum ~ 
        position_window * n_genes_Yunchen+position_window*expressed_leaf_sum +expressed_flower_sum + flower_female_bias_sum,
      #  n_genes_Yunchen,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      family = nbinom2,
      # ziformula=~1, 
      # control=glmmTMBControl(optimizer=optim,
      # optArgs=list(method="BFGS")),
        data = LG5)
    # data = LG5_distortion)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(female_cMMb_model5.4)
    anova(female_cMMb_model5.0, female_cMMb_model5.4)
    anova(female_cMMb_model5.1, female_cMMb_model5.4)
    anova(female_cMMb_model5.2, female_cMMb_model5.4)
    anova(female_cMMb_model5.3, female_cMMb_model5.4)
    
    
    simulationOutput <- simulateResiduals(fittedModel = female_cMMb_model5.4, plot = F, n = 1000)
    plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    plotResiduals(simulationOutput, LG5$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG5$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG5$sum_length_DNA_TEs_all)

    plotResiduals(simulationOutput, LG5$sum_length_RNA_TEs_all)

    plotResiduals(simulationOutput, LG5$female_distortion)
    plotResiduals(simulationOutput, LG5$male_distortion)
    plotResiduals(simulationOutput, LG5$pollen_biased_sum)
    plotResiduals(simulationOutput, LG5$pollen_tube_biased_sum)
    
    plotResiduals(simulationOutput,   (LG5$pollen_biased_sum/LG5$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG5$flower_female_bias_sum)
    plotResiduals(simulationOutput, LG5$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG5$expressed_flower_sum)
    
    
    plotResiduals(simulationOutput,   (LG5$flower_female_bias_sum/LG5$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG5$flower_female_bias_sum/LG5$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG5_distortion$sum_Length_Yunchen)
    plotResiduals(simulationOutput, LG5_distortion$n_genes_Yunchen)

    plotResiduals(simulationOutput, LG5_distortion$n_DNA_TEs_all)

    plotResiduals(simulationOutput, LG5_distortion$n_RNA_TEs_all)

    plotResiduals(simulationOutput, LG5_distortion$female_distortion)
    plotResiduals(simulationOutput, LG5_distortion$pollen_biased_sum)
    plotResiduals(simulationOutput, LG5_distortion$female_distortion_mean)
    
    plotResiduals(simulationOutput,   (LG5_distortion$pollen_biased_sum/LG5_distortion$n_genes_Yunchen))
    
    plotResiduals(simulationOutput, LG5_distortion$flower_female_bias_sum)
    plotResiduals(simulationOutput,   (LG5_distortion$flower_female_bias_sum/LG5_distortion$n_genes_Yunchen))
    plotResiduals(simulationOutput,   (LG5_distortion$flower_female_bias_sum/LG5_distortion$expressed_flower_sum))
    
    plotResiduals(simulationOutput, LG5_distortion$expressed_leaf_sum)
    plotResiduals(simulationOutput, LG5_distortion$expressed_flower_sum)
    
    
    ################ 4.2.3 Compare models for male- vs female-biased recombination  ################
    
    #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")

    
#    test <- subset(complete_windowed_dataset, complete_windowed_dataset$crossovers_diff == 0)
#    sum(test$crossovers_female_sum==0)


    
      complete_windowed_dataset <- complete_windowed_dataset %>%
      mutate(
        crossovers_binary = case_when(
        ((crossovers_male_sum > crossovers_female_sum) == TRUE) ~ 1,
        ((crossovers_male_sum < crossovers_female_sum) == TRUE) ~ 0))
    
    
    
    
    
        
    ################ 4.2.3.1 LG1 ################
    
    LG1<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG1")

 
    plot(LG1$position_window , LG1$crossovers_binary)
    plot(LG1$n_genes_Yunchen , LG1$crossovers_binary)
    plot(LG1$sum_Length_Yunchen , LG1$crossovers_binary)
    plot(LG1$female_distortion_mean, LG1$crossovers_binary)
    plot(LG1$N_DNA_TEs_all, LG1$crossovers_binary)
    plot(LG1$N_RNA_TEs_all, LG1$crossovers_binary)
    plot(LG1$pollen_biased_sum, LG1$crossovers_binary)
    plot(LG1$expressed_leaf_sum, LG1$crossovers_binary)
    plot(LG1$leaf_female_bias_sum, LG1$crossovers_binary)
    plot(jitter(LG1$flower_male_bias_sum), LG1$crossovers_binary)
    plot(jitter(LG1$flower_female_bias_sum), LG1$crossovers_binary)
    plot(LG1$pollen_tube_biased_sum, LG1$crossovers_binary)

    
    
    
    binary_crossover_1.0 <- glmmTMB(
      crossovers_binary ~ 
      #  position_window,
        n_genes_Yunchen,
      #family=gaussian,
      family=binomial,
#       ziformula=~1, 
      data = LG1)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_1.0)
    
    
    binary_crossover_1.1 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen*N_DNA_TEs_all,
        #,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG1)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_1.1)
    anova(binary_crossover_1.0,binary_crossover_1.1)
    
    
    binary_crossover_1.2 <- glmmTMB(
      crossovers_binary ~ 
      n_genes_Yunchen*N_RNA_TEs_all,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG1)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_1.2)
    anova(binary_crossover_1.0,binary_crossover_1.2)
    anova(binary_crossover_1.1,binary_crossover_1.2)
    
 
    
    binary_crossover_1.3 <- glmmTMB(
      crossovers_binary ~ 
      #  n_genes_Yunchen*N_RNA_TEs_all +leaf_male_bias_sum,
      n_genes_Yunchen*N_RNA_TEs_all+flower_female_bias_sum,
      #,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG1)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_1.3)
    anova(binary_crossover_1.0,binary_crossover_1.3)
    anova(binary_crossover_1.1,binary_crossover_1.3)
    anova(binary_crossover_1.2,binary_crossover_1.3)
    
       
    
    binary_crossover_1.4 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen*N_RNA_TEs_all+flower_female_bias_sum * female_distortion_mean,
      #,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG1)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_1.4)
    anova(binary_crossover_1.0,binary_crossover_1.4)
    anova(binary_crossover_1.1,binary_crossover_1.4)
    anova(binary_crossover_1.2,binary_crossover_1.4)
    anova(binary_crossover_1.3,binary_crossover_1.4)

    
    
    
    simulationOutput <- simulateResiduals(fittedModel = binary_crossover_1.4, plot = F, n = 1000)
    plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 

    
    
    ################ 4.2.3.2 LG2 ################
    
    LG2<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG2")
    
    
    plot(LG2$position_window , LG2$crossovers_binary)
    plot(LG2$n_genes_Yunchen , LG2$crossovers_binary)
    plot(LG2$sum_Length_Yunchen , LG2$crossovers_binary)
    plot(LG2$male_distortion_mean, LG2$crossovers_binary)
    plot(LG2$female_distortion_mean, LG2$crossovers_binary)
    plot(LG2$N_DNA_TEs_all, LG2$crossovers_binary)
    plot(LG2$N_RNA_TEs_all, LG2$crossovers_binary)
    plot(jitter(LG2$pollen_biased_sum), LG2$crossovers_binary)
    plot(jitter(LG2$expressed_leaf_sum), LG2$crossovers_binary)
    plot(jitter(LG2$leaf_female_bias_sum), LG2$crossovers_binary)
    plot(jitter(LG2$leaf_male_bias_sum), LG2$crossovers_binary)
    plot(jitter(LG2$flower_male_bias_sum), LG2$crossovers_binary)
    plot(jitter(LG2$flower_female_bias_sum), LG2$crossovers_binary)
    plot(jitter(LG2$pollen_tube_biased_sum), LG2$crossovers_binary)
    
    
    
    
    binary_crossover_2.0 <- glmmTMB(
      crossovers_binary ~ 
          position_window*n_genes_Yunchen,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG2)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_2.0)
    
    
    binary_crossover_2.1 <- glmmTMB(
      crossovers_binary ~ 
        position_window * N_DNA_TEs_all,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG2)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_2.1)
    anova(binary_crossover_2.0, binary_crossover_2.1)
    
    binary_crossover_2.2 <- glmmTMB(
      crossovers_binary ~ 
        position_window * N_DNA_TEs_all + N_DNA_TEs_all * expressed_leaf_sum,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG2)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_2.2)
    anova(binary_crossover_2.0, binary_crossover_2.2)
    anova(binary_crossover_2.1, binary_crossover_2.2)
    
 
    
    simulationOutput <- simulateResiduals(fittedModel = binary_crossover_2.2, plot = T, n = 1000)
    #plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 

    
    
    ################ 4.2.3.3 LG3 ################
    
    LG3<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG3")
    
    
    plot(LG3$position_window , LG3$crossovers_binary)
    plot(LG3$n_genes_Yunchen , LG3$crossovers_binary)
    plot(LG3$sum_Length_Yunchen , LG3$crossovers_binary)
    plot(LG3$male_distortion_mean, LG3$crossovers_binary)
    plot(LG3$female_distortion_mean, LG3$crossovers_binary)
    plot(LG3$N_DNA_TEs_all, LG3$crossovers_binary)
    plot(LG3$N_RNA_TEs_all, LG3$crossovers_binary)
    plot(jitter(LG3$pollen_biased_sum), LG3$crossovers_binary)
    plot(jitter(LG3$expressed_leaf_sum), LG3$crossovers_binary)
    plot(jitter(LG3$leaf_female_bias_sum), LG3$crossovers_binary)
    plot(jitter(LG3$leaf_male_bias_sum), LG3$crossovers_binary)
    plot(jitter(LG3$flower_male_bias_sum), LG3$crossovers_binary)
    plot(jitter(LG3$flower_female_bias_sum), LG3$crossovers_binary)
    plot(jitter(LG3$pollen_tube_biased_sum), LG3$crossovers_binary)
    
    
    
    
    binary_crossover_3.0 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG3)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_3.0)
    
    
    
    binary_crossover_3.1 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen * expressed_flower_sum,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG3)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_3.1)
    anova(binary_crossover_3.0,binary_crossover_3.1)
    
    binary_crossover_3.2 <- glmmTMB(
      crossovers_binary ~ 
      n_genes_Yunchen * flower_male_bias_sum,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG3)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_3.2)
    anova(binary_crossover_3.0,binary_crossover_3.2)
    anova(binary_crossover_3.1,binary_crossover_3.2)

    
    binary_crossover_3.3 <- glmmTMB(
      crossovers_binary ~ 
      n_genes_Yunchen * flower_male_bias_sum + 
      flower_male_bias_sum*male_distortion_mean , 
      family=binomial,
      #       ziformula=~1, 
      data = LG3)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_3.3)
    anova(binary_crossover_3.0, binary_crossover_3.3)
    anova(binary_crossover_3.1,binary_crossover_3.3)
    anova(binary_crossover_3.2,binary_crossover_3.3)
    
        
    
    simulationOutput <- simulateResiduals(fittedModel = binary_crossover_3.3, plot = T, n = 1000)
    #plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 

    
    
    ################ 4.2.3.4 LG4 ################
    
    LG4<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG4")
    
    
    plot(LG4$position_window , LG4$crossovers_binary)
    plot(LG4$n_genes_Yunchen , LG4$crossovers_binary)
    plot(LG4$sum_Length_Yunchen , LG4$crossovers_binary)
    plot(LG4$male_distortion_mean, LG4$crossovers_binary)
    plot(LG4$female_distortion_mean, LG4$crossovers_binary)
    plot(LG4$N_DNA_TEs_all, LG4$crossovers_binary)
    plot(LG4$N_RNA_TEs_all, LG4$crossovers_binary)
    plot(jitter(LG4$pollen_biased_sum), LG4$crossovers_binary)
    plot(jitter(LG4$expressed_leaf_sum), LG4$crossovers_binary)
    plot(jitter(LG4$leaf_female_bias_sum), LG4$crossovers_binary)
    plot(jitter(LG4$leaf_male_bias_sum), LG4$crossovers_binary)
    plot(jitter(LG4$flower_male_bias_sum), LG4$crossovers_binary)
    plot(jitter(LG4$flower_female_bias_sum), LG4$crossovers_binary)
    plot(jitter(LG4$pollen_tube_biased_sum), LG4$crossovers_binary)
    
    
    binary_crossover_4.0 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen,
       #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG4)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_4.0)
    
    binary_crossover_4.1 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen * expressed_leaf_sum,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG4)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_4.1)
    anova(binary_crossover_4.0 , binary_crossover_4.1)
    
    
    binary_crossover_4.2 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen * expressed_leaf_sum + flower_female_bias_sum,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG4)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_4.2)
    anova(binary_crossover_4.0 , binary_crossover_4.2)
    anova(binary_crossover_4.1 , binary_crossover_4.2)
    
    
    binary_crossover_4.3 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen * expressed_flower_sum + n_genes_Yunchen* pollen_biased_sum,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG4)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_4.3)
    anova(binary_crossover_4.0 , binary_crossover_4.3)
    anova(binary_crossover_4.1 , binary_crossover_4.3)
    anova(binary_crossover_4.2 , binary_crossover_4.3)
    

    
    
    simulationOutput <- simulateResiduals(fittedModel = binary_crossover_4.3, plot = T, n = 1000)
    #plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 

    
    ################ 4.2.3.5 LG5 ################
    
    LG5<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG5")
    
    
    plot(LG5$position_window , LG5$crossovers_binary)
    plot(LG5$n_genes_Yunchen , LG5$crossovers_binary)
    plot(LG5$sum_Length_Yunchen , LG5$crossovers_binary)
    plot(LG5$male_distortion_mean, LG5$crossovers_binary)
    plot(LG5$female_distortion_mean, LG5$crossovers_binary)
    plot(LG5$N_DNA_TEs_all, LG5$crossovers_binary)
    plot(LG5$N_RNA_TEs_all, LG5$crossovers_binary)
    plot(jitter(LG5$pollen_biased_sum), LG5$crossovers_binary)
    plot(jitter(LG5$expressed_leaf_sum), LG5$crossovers_binary)
    plot(jitter(LG5$leaf_female_bias_sum), LG5$crossovers_binary)
    plot(jitter(LG5$leaf_male_bias_sum), LG5$crossovers_binary)
    plot(jitter(LG5$flower_male_bias_sum), LG5$crossovers_binary)
    plot(jitter(LG5$flower_female_bias_sum), LG5$crossovers_binary)
    plot(jitter(LG5$pollen_tube_biased_sum), LG5$crossovers_binary)
    
    
    binary_crossover_5.0 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen*position_window,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG5)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_5.0)
    
    binary_crossover_5.1 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen*position_window + pollen_tube_biased_sum,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG5)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_5.1)
    anova(binary_crossover_5.0,binary_crossover_5.1)
    
    binary_crossover_5.2 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen*position_window + pollen_tube_biased_sum + position_window*female_distortion_mean,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG5)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_5.2)
    anova(binary_crossover_5.0,binary_crossover_5.2)
    anova(binary_crossover_5.1,binary_crossover_5.2)
    
    binary_crossover_5.3 <- glmmTMB(
      crossovers_binary ~ 
        n_genes_Yunchen*position_window + position_window*female_distortion_mean,
      #family=gaussian,
      family=binomial,
      #       ziformula=~1, 
      data = LG5)
    #negative binomial slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(binary_crossover_5.3)
    anova(binary_crossover_5.0,binary_crossover_5.3)
    anova(binary_crossover_5.1,binary_crossover_5.3)
    anova(binary_crossover_5.2,binary_crossover_5.3)
    
    
    
        
    simulationOutput <- simulateResiduals(fittedModel = binary_crossover_5.3, plot = T, n = 1000)
    #plot(simulationOutput)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 

    
    
    ################ 4.2.4 Compare models for recombination rate difference ################
    #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
    
    complete_windowed_dataset$abs_crossovers_diff<-abs(complete_windowed_dataset$crossovers_male_sum-complete_windowed_dataset$crossovers_female_sum)
    
    
    ################ 4.2.4.1 LG1 ################
    
    
    LG1<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG1")
    LG1_distortion<-subset(LG1, !is.na(female_distortion_mean))
    
    
    plot(LG1$position_window , LG1$abs_crossovers_diff)
    plot(LG1$n_genes_Yunchen , LG1$abs_crossovers_diff)
    plot(LG1$sum_Length_Yunchen , LG1$abs_crossovers_diff)
    plot(LG1$female_distortion_mean, LG1$abs_crossovers_diff)
    plot(LG1$N_DNA_TEs_all, LG1$abs_crossovers_diff)
    plot(LG1$N_RNA_TEs_all, LG1$abs_crossovers_diff)
    plot(LG1$pollen_biased_sum, LG1$abs_crossovers_diff)
    plot(LG1$expressed_leaf_sum, LG1$abs_crossovers_diff)
    plot(LG1$leaf_female_bias_sum, LG1$abs_crossovers_diff)
    plot(jitter(LG1$flower_male_bias_sum), LG1$abs_crossovers_diff)
    plot(jitter(LG1$flower_female_bias_sum), LG1$abs_crossovers_diff)
    plot(LG1$pollen_tube_biased_sum, LG1$abs_crossovers_diff)
    
    
    
    
    crossover_diff_1.0 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen,
       family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
     #        ziformula=~1, 
   #  data = LG1)
    data = LG1_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(crossover_diff_1.0)
    
    
    
    crossover_diff_1.1 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen*N_RNA_TEs_all,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #        ziformula=~1, 
    #  data = LG1)
      data = LG1_distortion)
    #run on both full and distortion-only dataset
    summary(crossover_diff_1.1)
    anova(crossover_diff_1.0, crossover_diff_1.1)
    
    crossover_diff_1.2 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen*pollen_tube_biased_sum+n_genes_Yunchen*N_RNA_TEs_all+N_RNA_TEs_all*pollen_tube_biased_sum,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #        ziformula=~1, 
     #   data = LG1)
      data = LG1_distortion)
    #run on both full and distortion-only dataset
    summary(crossover_diff_1.2)
    anova(crossover_diff_1.0, crossover_diff_1.2)
    anova(crossover_diff_1.1, crossover_diff_1.2)

    
    crossover_diff_1.3 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen*pollen_tube_biased_sum+n_genes_Yunchen*N_RNA_TEs_all+N_RNA_TEs_all*pollen_tube_biased_sum+male_distortion_mean+female_distortion_mean,        
        #        n_genes_Yunchen*pollen_tube_biased_sum+n_genes_Yunchen*N_RNA_TEs_all+N_RNA_TEs_all*pollen_tube_biased_sum+male_distortion_mean,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #        ziformula=~1, 
         data = LG1)
      #data = LG1_distortion)
    #run on both full and distortion-only dataset
    summary(crossover_diff_1.3)
    anova(crossover_diff_1.0, crossover_diff_1.2)
    anova(crossover_diff_1.1, crossover_diff_1.2)
    
    
    
        
    simulationOutput <- simulateResiduals(fittedModel = crossover_diff_1.2, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG1
    testdata<-LG1_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$sum_length_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$sum_length_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)
    
 
    ################ 4.2.4.2 LG2 ################
    
    
    LG2<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG2")
    LG2_distortion<-subset(LG2, !is.na(female_distortion_mean))
    
    hist(LG2$abs_crossovers_diff)
    
    plot(LG2$position_window , LG2$abs_crossovers_diff)
    plot(LG2$n_genes_Yunchen , LG2$abs_crossovers_diff)
    plot(LG2$sum_Length_Yunchen , LG2$abs_crossovers_diff)
    plot(LG2$male_distortion_mean, LG2$abs_crossovers_diff)
    plot(LG2$female_distortion_mean, LG2$abs_crossovers_diff)
    plot(LG2$N_DNA_TEs_all, LG2$abs_crossovers_diff)
    plot(LG2$N_RNA_TEs_all, LG2$abs_crossovers_diff)
    plot(LG2$pollen_biased_sum, LG2$abs_crossovers_diff)
    plot(LG2$expressed_leaf_sum, LG2$abs_crossovers_diff)
    plot(LG2$leaf_female_bias_sum, LG2$abs_crossovers_diff)
    plot(jitter(LG2$flower_male_bias_sum), LG2$abs_crossovers_diff)
    plot(jitter(LG2$flower_female_bias_sum), LG2$abs_crossovers_diff)
    plot(LG2$pollen_biased_sum, LG2$abs_crossovers_diff)
    plot(LG2$pollen_tube_biased_sum, LG2$abs_crossovers_diff)
    
    
    
    
    crossover_diff_2.0 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen*position_window,
      family = nbinom2,
    #    family = poisson(link = "log"), 
     #  family = tweedie(link = "log"), 
      #  ziformula=~1, 
        data = LG2)
    #  data = LG2_distortion)
     but overdispersed


    #run on both full and distortion-only dataset
    summary(crossover_diff_2.0)
    

    crossover_diff_2.1 <- glmmTMB(
      abs_crossovers_diff ~ 
      n_genes_Yunchen*position_window+position_window*flower_female_bias_sum,
      family = nbinom2,
      #    family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      #  ziformula=~1, 
      data = LG2)
    #  data = LG2_distortion)

    #run on both full and distortion-only dataset
    summary(crossover_diff_2.1)
    anova(crossover_diff_2.0, crossover_diff_2.1)

    
    
    crossover_diff_2.2 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen*position_window+flower_female_bias_sum,
      family = nbinom2,
      #    family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      #  ziformula=~1, 
    #  data = LG2)
      data = LG2_distortion)
    


    #run on both full and distortion-only dataset
    summary(crossover_diff_2.2)
    anova(crossover_diff_2.0, crossover_diff_2.2)
    anova(crossover_diff_2.1, crossover_diff_2.2)
  
    simulationOutput <- simulateResiduals(fittedModel = crossover_diff_2.2, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG2
    testdata<-LG2_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$sum_length_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$sum_length_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)
    
    ################ 4.2.4.3 LG3 ################
    
    
    LG3<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG3")
    LG3_distortion<-subset(LG3, !is.na(female_distortion_mean))
    
    hist(LG3$abs_crossovers_diff)
    
    plot(LG3$position_window , LG3$abs_crossovers_diff)
    plot(LG3$n_genes_Yunchen , LG3$abs_crossovers_diff)
    plot(LG3$sum_Length_Yunchen , LG3$abs_crossovers_diff)
    plot(LG3$male_distortion_mean, LG3$abs_crossovers_diff)
    plot(LG3$female_distortion_mean, LG3$abs_crossovers_diff)
    plot(LG3$N_DNA_TEs_all, LG3$abs_crossovers_diff)
    plot(LG3$N_RNA_TEs_all, LG3$abs_crossovers_diff)
    plot(LG3$pollen_biased_sum, LG3$abs_crossovers_diff)
    plot(LG3$expressed_leaf_sum, LG3$abs_crossovers_diff)
    plot(LG3$leaf_female_bias_sum, LG3$abs_crossovers_diff)
    plot(jitter(LG3$flower_male_bias_sum), LG3$abs_crossovers_diff)
    plot(jitter(LG3$flower_female_bias_sum), LG3$abs_crossovers_diff)
    plot(LG3$pollen_biased_sum, LG3$abs_crossovers_diff)
    plot(LG3$pollen_tube_biased_sum, LG3$abs_crossovers_diff)
    
    
    
    
    crossover_diff_3.0 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen+position_window,
      family = nbinom2,
      #    family = poisson(link = "log"), 
       # family = tweedie(link = "log"), 
      #  ziformula=~1, 
      data = LG3)
    #  data = LG3_distortion)
    #negative binomial without zero-inflation -  underdispersed


    #run on both full and distortion-only dataset
    summary(crossover_diff_3.0)
    
    
    crossover_diff_3.0 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen+position_window,
      family = nbinom2,
      #    family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #  ziformula=~1, 
   #   data = LG3)
      data = LG3_distortion)
    #negative binomial without zero-inflation -  underdispersed


    #run on both full and distortion-only dataset
    summary(crossover_diff_3.0)
    
    
    crossover_diff_3.1 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen*N_RNA_TEs_all+position_window ,
      family = nbinom2,
      #    family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #  ziformula=~1, 
    #  data = LG3)
      data = LG3_distortion)
    #negative binomial without zero-inflation -  underdispersed


    #run on both full and distortion-only dataset
    summary(crossover_diff_3.1)
    anova(crossover_diff_3.0, crossover_diff_3.1)
    
    
    crossover_diff_3.2 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen*N_RNA_TEs_all+n_genes_Yunchen*male_distortion_mean,
      family = nbinom2,
      #    family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #  ziformula=~1, 
   #   data = LG3)
      data = LG3_distortion)
    #negative binomial without zero-inflation -  underdispersed


    #run on both full and distortion-only dataset
    summary(crossover_diff_3.2)
    anova(crossover_diff_3.0, crossover_diff_3.2)
    anova(crossover_diff_3.1, crossover_diff_3.2)
    
    
    simulationOutput <- simulateResiduals(fittedModel = crossover_diff_3.1, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG3
    testdata<-LG3_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$sum_length_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$sum_length_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)
    
    ################ 4.2.4.4 LG4 ################
    
    
    LG4<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG4")
    LG4_distortion<-subset(LG4, !is.na(female_distortion_mean))
    
    hist(LG4$abs_crossovers_diff)
    
    plot(LG4$position_window , LG4$abs_crossovers_diff)
    plot(LG4$n_genes_Yunchen , LG4$abs_crossovers_diff)
    plot(LG4$sum_Length_Yunchen , LG4$abs_crossovers_diff)
    plot(LG4$male_distortion_mean, LG4$abs_crossovers_diff)
    plot(LG4$female_distortion_mean, LG4$abs_crossovers_diff)
    plot(LG4$N_DNA_TEs_all, LG4$abs_crossovers_diff)
    plot(LG4$N_RNA_TEs_all, LG4$abs_crossovers_diff)
    plot(LG4$pollen_biased_sum, LG4$abs_crossovers_diff)
    plot(LG4$expressed_leaf_sum, LG4$abs_crossovers_diff)
    plot(LG4$leaf_female_bias_sum, LG4$abs_crossovers_diff)
    plot(jitter(LG4$flower_male_bias_sum), LG4$abs_crossovers_diff)
    plot(jitter(LG4$flower_female_bias_sum), LG4$abs_crossovers_diff)
    plot(LG4$pollen_biased_sum, LG4$abs_crossovers_diff)
    plot(LG4$pollen_tube_biased_sum, LG4$abs_crossovers_diff)
    
    
    
    
    crossover_diff_4.0 <- glmmTMB(
      abs_crossovers_diff ~ 
        n_genes_Yunchen+position_window,
      family = nbinom2,
      #    family = poisson(link = "log"), 
     #  family = tweedie(link = "log"), 
      #  ziformula=~1, 
   #   data = LG4)
      data = LG4_distortion)
    #tweedie fits better
    #run on both full and distortion-only dataset
    summary(crossover_diff_4.0)   
    
    crossover_diff_4.1 <- glmmTMB(
      abs_crossovers_diff ~ 
          n_genes_Yunchen+position_window*flower_female_bias_sum,
  family = nbinom2,
      #    family = poisson(link = "log"), 
      #family = tweedie(link = "log"), 
      #  ziformula=~1, 
    #  data = LG4)
      data = LG4_distortion)
    #run on both full and distortion-only dataset
    summary(crossover_diff_4.1)  
    anova(crossover_diff_4.0, crossover_diff_4.1)


    crossover_diff_4.2 <- glmmTMB(
      abs_crossovers_diff ~ 
         position_window*female_distortion_mean,
  family = nbinom2,
      #    family = poisson(link = "log"), 
      #family = tweedie(link = "log"), 
      #  ziformula=~1, 
    #  data = LG4)
      data = LG4_distortion)
    summary(crossover_diff_4.2)  
    anova(crossover_diff_4.0, crossover_diff_4.2)
    anova(crossover_diff_4.1, crossover_diff_4.2)
    
    
            
    simulationOutput <- simulateResiduals(fittedModel = crossover_diff_4.2, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG4
    testdata<-LG4_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$sum_length_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$sum_length_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)
    
    
    ################ 4.2.4.5 LG5 ################
    
    
    LG5<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG5")
    LG5_distortion<-subset(LG5, !is.na(female_distortion_mean))
    
    hist(LG5$abs_crossovers_diff)
    
    plot(LG5$position_window , LG5$abs_crossovers_diff)
    plot(LG5$n_genes_Yunchen , LG5$abs_crossovers_diff)
    plot(LG5$sum_Length_Yunchen , LG5$abs_crossovers_diff)
    plot(LG5$male_distortion_mean, LG5$abs_crossovers_diff)
    plot(LG5$female_distortion_mean, LG5$abs_crossovers_diff)
    plot(LG5$N_DNA_TEs_all, LG5$abs_crossovers_diff)
    plot(LG5$N_RNA_TEs_all, LG5$abs_crossovers_diff)
    plot(LG5$pollen_biased_sum, LG5$abs_crossovers_diff)
    plot(LG5$expressed_leaf_sum, LG5$abs_crossovers_diff)
    plot(LG5$leaf_female_bias_sum, LG5$abs_crossovers_diff)
    plot(jitter(LG5$flower_male_bias_sum), LG5$abs_crossovers_diff)
    plot(jitter(LG5$flower_female_bias_sum), LG5$abs_crossovers_diff)
    plot(LG5$pollen_biased_sum, LG5$abs_crossovers_diff)
    plot(LG5$pollen_tube_biased_sum, LG5$abs_crossovers_diff)
    
    
    
    
    crossover_diff_5.0 <- glmmTMB(
      abs_crossovers_diff ~ 
        position_window,
      family = nbinom2,
    #      family = poisson(link = "log"), 
     #   family = tweedie(link = "log"), 
      #  ziformula=~1, 
  #      data = LG5)
      data = LG5_distortion)
    #run on both full and distortion-only dataset
    summary(crossover_diff_5.0)   

    
    crossover_diff_5.1 <- glmmTMB(
      abs_crossovers_diff ~ 
        position_window*male_distortion_mean,
      family = nbinom2,
      #      family = poisson(link = "log"), 
      #   family = tweedie(link = "log"), 
      #  ziformula=~1, 
    #  data = LG5)
          data = LG5_distortion)
    #run on both full and distortion-only dataset
    summary(crossover_diff_5.1)   
    anova(crossover_diff_5.0,crossover_diff_5.1)   
    
    crossover_diff_5.2 <- glmmTMB(
      abs_crossovers_diff ~ 
        position_window*male_distortion_mean+flower_female_bias_sum,
      family = nbinom2,
      #      family = poisson(link = "log"), 
      #   family = tweedie(link = "log"), 
      #  ziformula=~1, 
     #   data = LG5)
      data = LG5_distortion)
    #run on both full and distortion-only dataset
    summary(crossover_diff_5.2)   
    anova(crossover_diff_5.0,crossover_diff_5.2)   
    anova(crossover_diff_5.1,crossover_diff_5.2)   

    
    simulationOutput <- simulateResiduals(fittedModel = crossover_diff_5.2, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG5
    testdata<-LG5_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$sum_length_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$sum_length_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)
    
     
    
    ################ 4.2.5 Compare models for total recombination  ################
    #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
    
    complete_windowed_dataset$total_crossovers_sum<-(complete_windowed_dataset$crossovers_sex_averaged_sum)*2

    ################ 4.2.5.1 LG1 ################
    
    
    LG1<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG1")
    LG1_distortion<-subset(LG1, !is.na(female_distortion_mean))
    
    plot(LG1$position_window , LG1$total_crossovers_sum)
    plot(LG1$n_genes_Yunchen , LG1$total_crossovers_sum)
    plot(LG1$sum_Length_Yunchen , LG1$total_crossovers_sum)
    plot(LG1$female_distortion_mean, LG1$total_crossovers_sum)
    plot(LG1$N_DNA_TEs_all, LG1$total_crossovers_sum)
    plot(LG1$N_RNA_TEs_all, LG1$total_crossovers_sum)
    plot(LG1$pollen_biased_sum, LG1$total_crossovers_sum)
    plot(LG1$expressed_leaf_sum, LG1$total_crossovers_sum)
    plot(LG1$leaf_female_bias_sum, LG1$total_crossovers_sum)
    plot(jitter(LG1$flower_male_bias_sum), LG1$total_crossovers_sum)
    plot(jitter(LG1$flower_female_bias_sum), LG1$total_crossovers_sum)
    plot(LG1$pollen_tube_biased_sum, LG1$total_crossovers_sum)
    

    sex_avg_1.0 <- glmmTMB(
      total_crossovers_sum ~ 
      n_genes_Yunchen ,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
       #       ziformula=~1, 
    #    data = LG1)
      data = LG1_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_1.0)
       

    
    sex_avg_1.1 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen * N_RNA_TEs_all,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
   #  data = LG1)
   data = LG1_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_1.1)
    anova(sex_avg_1.0,sex_avg_1.1)
    
    
    sex_avg_1.2 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen  + N_RNA_TEs_all*pollen_tube_biased_sum,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
   #   data = LG1)
    data = LG1_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_1.2)
    anova(sex_avg_1.0,sex_avg_1.2)
    anova(sex_avg_1.1,sex_avg_1.2)
    
    
    sex_avg_1.3 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen*pollen_tube_biased_sum,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
    #     data = LG1)
     data = LG1_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_1.3)
    anova(sex_avg_1.0,sex_avg_1.3)
    anova(sex_avg_1.1,sex_avg_1.3)
    anova(sex_avg_1.2,sex_avg_1.3)
    
    
    sex_avg_1.4 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen + pollen_tube_biased_sum +n_genes_Yunchen*male_distortion_mean,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
      #   data = LG1)
      data = LG1_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_1.4)
    anova(sex_avg_1.0,sex_avg_1.4)
    anova(sex_avg_1.1,sex_avg_1.4)
    anova(sex_avg_1.2,sex_avg_1.4)
    anova(sex_avg_1.3,sex_avg_1.4)
    
    
    sex_avg_1.5 <- glmmTMB(
      total_crossovers_sum ~ 
      pollen_tube_biased_sum +n_genes_Yunchen*male_distortion_mean + female_distortion_mean,
    family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
      #   data = LG1)
      data = LG1_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_1.5)
    anova(sex_avg_1.0,sex_avg_1.5)
    anova(sex_avg_1.1,sex_avg_1.5)
    anova(sex_avg_1.2,sex_avg_1.5)
    anova(sex_avg_1.3,sex_avg_1.5)
    anova(sex_avg_1.4,sex_avg_1.5)
    
    
    
    
    simulationOutput <- simulateResiduals(fittedModel = sex_avg_1.5, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG1
    testdata<-LG1_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$N_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$N_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)
    
    ################ 4.2.5.2 LG2 ################
    
    
    LG2<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG2")
    LG2_distortion<-subset(LG2, !is.na(female_distortion_mean))
    
    plot(LG2$position_window , LG2$total_crossovers_sum)
    plot(LG2$n_genes_Yunchen , LG2$total_crossovers_sum)
    plot(LG2$sum_Length_Yunchen , LG2$total_crossovers_sum)
    plot(LG2$female_distortion_mean, LG2$total_crossovers_sum)
    plot(LG2$N_DNA_TEs_all, LG2$total_crossovers_sum)
    plot(LG2$N_RNA_TEs_all, LG2$total_crossovers_sum)
    plot(LG2$pollen_biased_sum, LG2$total_crossovers_sum)
    plot(LG2$expressed_leaf_sum, LG2$total_crossovers_sum)
    plot(LG2$leaf_female_bias_sum, LG2$total_crossovers_sum)
    plot(jitter(LG2$flower_male_bias_sum), LG2$total_crossovers_sum)
    plot(jitter(LG2$flower_female_bias_sum), LG2$total_crossovers_sum)
    plot(jitter(LG2$pollen_tube_biased_sum), LG2$total_crossovers_sum)
    
    
    sex_avg_2.0 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen *position_window,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
     #     data = LG2)
      data = LG2_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_2.0)
    
    
    sex_avg_2.1 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen *position_window + position_window *flower_female_bias_sum,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
   #   data = LG2)
    data = LG2_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_2.1)
    anova(sex_avg_2.0,sex_avg_2.1)
    
        
    sex_avg_2.2 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen *position_window +male_distortion_mean,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
      #   data = LG2)
      data = LG2_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_2.2)
    anova(sex_avg_2.0,sex_avg_2.2)
    anova(sex_avg_2.1,sex_avg_2.2)
    
    
    
        
    simulationOutput <- simulateResiduals(fittedModel = sex_avg_2.2, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG2
    testdata<-LG2_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$N_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$N_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)

    
    
    
    ################ 4.2.5.3 LG3 ################
    
    
    LG3<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG3")
    LG3_distortion<-subset(LG3, !is.na(female_distortion_mean))
    
    plot(LG3$position_window , LG3$total_crossovers_sum)
    plot(LG3$n_genes_Yunchen , LG3$total_crossovers_sum)
    plot(LG3$sum_Length_Yunchen , LG3$total_crossovers_sum)
    plot(LG3$female_distortion_mean, LG3$total_crossovers_sum)
    plot(LG3$N_DNA_TEs_all, LG3$total_crossovers_sum)
    plot(LG3$N_RNA_TEs_all, LG3$total_crossovers_sum)
    plot(LG3$pollen_biased_sum, LG3$total_crossovers_sum)
    plot(LG3$expressed_leaf_sum, LG3$total_crossovers_sum)
    plot(LG3$leaf_female_bias_sum, LG3$total_crossovers_sum)
    plot(jitter(LG3$flower_male_bias_sum), LG3$total_crossovers_sum)
    plot(jitter(LG3$flower_female_bias_sum), LG3$total_crossovers_sum)
    plot(jitter(LG3$pollen_tube_biased_sum), LG3$total_crossovers_sum)
    
    
    sex_avg_3.0 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen *position_window,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
      #     data = LG3)
     data = LG3_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_3.0)


    
    sex_avg_3.1 <- glmmTMB(
      total_crossovers_sum ~ 
       n_genes_Yunchen *position_window+n_genes_Yunchen *pollen_biased_sum,
        family = nbinom2,
     #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
     #       data = LG3)
      data = LG3_distortion)
    #negative binomial, slightly better without zero-inflation


    #run on both full and distortion-only dataset
    summary(sex_avg_3.1)
    anova(sex_avg_3.0,sex_avg_3.1)
    
   
    simulationOutput <- simulateResiduals(fittedModel = sex_avg_3.1, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG3
    testdata<-LG3_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$N_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$N_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)
    
    
    ################ 4.2.5.4 LG4 ################
    
    
    LG4<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG4")
    LG4_distortion<-subset(LG4, !is.na(female_distortion_mean))
    
    plot(LG4$position_window , LG4$total_crossovers_sum)
    plot(LG4$n_genes_Yunchen , LG4$total_crossovers_sum)
    plot(LG4$sum_Length_Yunchen , LG4$total_crossovers_sum)
    plot(LG4$female_distortion_mean, LG4$total_crossovers_sum)
    plot(LG4$N_DNA_TEs_all, LG4$total_crossovers_sum)
    plot(LG4$N_RNA_TEs_all, LG4$total_crossovers_sum)
    plot(LG4$pollen_biased_sum, LG4$total_crossovers_sum)
    plot(LG4$expressed_leaf_sum, LG4$total_crossovers_sum)
    plot(LG4$leaf_female_bias_sum, LG4$total_crossovers_sum)
    plot(jitter(LG4$flower_male_bias_sum), LG4$total_crossovers_sum)
    plot(jitter(LG4$flower_female_bias_sum), LG4$total_crossovers_sum)
    plot(jitter(LG4$pollen_tube_biased_sum), LG4$total_crossovers_sum)
    
    
    
    sex_avg_4.0 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen +position_window,
     # family = nbinom2,
      #  family = poisson(link = "log"), 
       family = tweedie(link = "log"), 
      #       ziformula=~1, 
           data = LG4)
      #data = LG4_distortion)
    #tweedie fits better than negative binomial 
    summary(sex_avg_4.0)
            
    
    sex_avg_4.0 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen +position_window,
       family = nbinom2,
      #  family = poisson(link = "log"), 
    #  family = tweedie(link = "log"), 
      #       ziformula=~1, 
      #data = LG4)
    data = LG4_distortion)
    #negative binomial fits better with distortion dataset
    summary(sex_avg_4.0)
    
    
    sex_avg_4.1 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen +position_window+position_window*female_distortion_mean,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      #  family = tweedie(link = "log"), 
      #       ziformula=~1, 
      #data = LG4)
      data = LG4_distortion)
    #run on both full and distortion-only dataset
    summary(sex_avg_4.1)
    anova(sex_avg_4.0,sex_avg_4.1)
     
    sex_avg_4.2 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen +position_window + flower_female_bias_sum,
      # family = nbinom2,
      #  family = poisson(link = "log"), 
      family = tweedie(link = "log"), 
      #       ziformula=~1, 
      data = LG4)
    #data = LG4_distortion)
    #tweedie fits better than negative binomial 
    #run on both full and distortion-only dataset
    summary(sex_avg_4.2)
    anova(sex_avg_4.0,sex_avg_4.2)
    
    
    sex_avg_4.3 <- glmmTMB(
      total_crossovers_sum ~ 
        n_genes_Yunchen +position_window + position_window*flower_female_bias_sum,
      # family = nbinom2,
      #  family = poisson(link = "log"), 
      family = tweedie(link = "log"), 
      #       ziformula=~1, 
      data = LG4)
    #data = LG4_distortion)
    #tweedie fits better than negative binomial 
    #run on both full and distortion-only dataset
    summary(sex_avg_4.3)
    anova(sex_avg_4.0,sex_avg_4.3)
    anova(sex_avg_4.2,sex_avg_4.3)
    
 
    
    simulationOutput <- simulateResiduals(fittedModel = sex_avg_4.3, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG4
    testdata<-LG4_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$N_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$N_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)
    
   
    
    
    
    ################ 4.2.5.4 LG5 ################
    
    
    LG5<-subset(complete_windowed_dataset,complete_windowed_dataset$LG=="LG5")
    LG5_distortion<-subset(LG5, !is.na(female_distortion_mean))
    
    plot(LG5$position_window , LG5$total_crossovers_sum)
    plot(LG5$n_genes_Yunchen , LG5$total_crossovers_sum)
    plot(LG5$sum_Length_Yunchen , LG5$total_crossovers_sum)
    plot(LG5$female_distortion_mean, LG5$total_crossovers_sum)
    plot(LG5$N_DNA_TEs_all, LG5$total_crossovers_sum)
    plot(LG5$N_RNA_TEs_all, LG5$total_crossovers_sum)
    plot(LG5$pollen_biased_sum, LG5$total_crossovers_sum)
    plot(LG5$expressed_leaf_sum, LG5$total_crossovers_sum)
    plot(LG5$leaf_female_bias_sum, LG5$total_crossovers_sum)
    plot(jitter(LG5$flower_male_bias_sum), LG5$total_crossovers_sum)
    plot(jitter(LG5$flower_female_bias_sum), LG5$total_crossovers_sum)
    plot(jitter(LG5$pollen_tube_biased_sum), LG5$total_crossovers_sum)
    
    
    
    
    
    
    
    sex_avg_5.0 <- glmmTMB(
      total_crossovers_sum ~ 
        position_window,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
    #  data = LG5)
    data = LG5_distortion)
    #tweedie fits better than negative binomial 
    summary(sex_avg_5.0)
    
    
    sex_avg_5.1 <- glmmTMB(
      total_crossovers_sum ~ 
        position_window*male_distortion_mean,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
         #data = LG5)
      data = LG5_distortion)
    #tweedie fits better than negative binomial 
    summary(sex_avg_5.1)
    anova(sex_avg_5.0,sex_avg_5.1)
  
    sex_avg_5.2 <- glmmTMB(
      total_crossovers_sum ~ 
        position_window*female_distortion_mean + position_window*male_distortion_mean,
      family = nbinom2,
      #  family = poisson(link = "log"), 
      # family = tweedie(link = "log"), 
      #       ziformula=~1, 
      #data = LG5)
      data = LG5_distortion)
    #tweedie fits better than negative binomial 
    summary(sex_avg_5.2)
    anova(sex_avg_5.0,sex_avg_5.2)
    anova(sex_avg_5.1,sex_avg_5.2)
    
      
    simulationOutput <- simulateResiduals(fittedModel = sex_avg_5.2, plot = T, n = 1000)




    testQuantiles(simulationOutput)
    testDispersion(simulationOutput)

    testZeroInflation(simulationOutput) 


    testdata<-LG5
    testdata<-LG5_distortion
    plotResiduals(simulationOutput, testdata$sum_Length_Yunchen)
    plotResiduals(simulationOutput, testdata$n_genes_Yunchen)
    plotResiduals(simulationOutput, testdata$N_DNA_TEs_all)
    plotResiduals(simulationOutput, testdata$N_RNA_TEs_all)
    plotResiduals(simulationOutput, testdata$male_distortion_mean)
    plotResiduals(simulationOutput, testdata$female_distortion_mean)
    plotResiduals(simulationOutput, testdata$expressed_leaf_sum)
    plotResiduals(simulationOutput, testdata$expressed_flower_sum)
    plotResiduals(simulationOutput, testdata$leaf_male_bias_sum)
    plotResiduals(simulationOutput, testdata$leaf_female_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_male_bias_sum)
    plotResiduals(simulationOutput, testdata$flower_female_bias_sum)
    plotResiduals(simulationOutput, testdata$pollen_biased_sum)
    plotResiduals(simulationOutput, testdata$pollen_tube_biased_sum)
    
    
    
    
    
    
################ 5. Create figures ################
     
    # fullmap<- read.csv("Windowed_analyses/fullmap_1.4_9-10-2021.csv")
    #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
    
    plotting_map<-fullmap %>% dplyr::select(LG, bp_position, cm_male, cm_female, cm_sex_averaged, crossovers_male, crossovers_female, male_distortion,female_distortion)
    
    plotting_map$mb_position<-plotting_map$bp_position/1000000

    #Perfect sizing Fig1: https://stackoverflow.com/questions/49110877/how-to-adjust-facet-size-manually
    #Panel tweaking: https://cran.r-project.org/web/packages/egg/vignettes/Overview.html

    ################ 5.1 Figure 1 ################
    
    ################ 5.1.1 Fig 1A Marey maps ################
    reshaped_map<-plotting_map%>%pivot_longer(!c(LG,bp_position), names_to = "chr", values_to = "cM")
    
    reshaped_map<-fullmap%>%pivot_longer(c(cm_male, cm_female), names_to = "sex", values_to = "cM")
    
    reshaped_map$LG[reshaped_map$LG == "LG4"] <- "XY"
    reshaped_map$LG[reshaped_map$LG == "LG1"] <- "A1"
    reshaped_map$LG[reshaped_map$LG == "LG3"] <- "A2"
    reshaped_map$LG[reshaped_map$LG == "LG5"] <- "A3"
    reshaped_map$LG[reshaped_map$LG == "LG2"] <- "A4"
    reshaped_map$LG = factor(reshaped_map$LG, levels=c('A1','A2','A3','A4','XY'))
    
   # sampleNo = 10000
    #reshaped_map <- reshaped_map[sample(nrow(reshaped_map), sampleNo),]

    marey_plot_m_f <-
      ggplot(reshaped_map) + 
      geom_point(alpha=0.5,aes(x=bp_position,y=cM,color=sex)) + 
      #  geom_line(aes(x=mb,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
      facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
      theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
      labs(x="",y="Genetic distance (cM)",color="Sex") +
      theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      scale_color_manual(values=c( 
        'darkgreen',
        'purple' #Blue #X
      )) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none") 
     #theme(axis.title.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))
    

    marey_A1<-ggplot(subset(reshaped_map, reshaped_map$LG=="A1")) + 
      geom_line(aes(x=bp_position,y=cM,color=sex), lineend="round", linejoin="round") +
      theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
#      geom_point(alpha=0.5, size=1, aes(x=bp_position,y=cM,color=sex)) +
#      theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
      #labs(x="",y="Genetic distance (cM)",color="Sex") +
      theme(strip.background =element_rect(fill="white")) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      scale_color_manual(values=c( 
        'darkgreen',
        'purple' #Blue #X
      )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none") 
    
    
    marey_A2<-ggplot(subset(reshaped_map, reshaped_map$LG=="A2")) + 
      geom_line(aes(x=bp_position,y=cM,color=sex), lineend="round", linejoin="round") +
      theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
#      geom_point(alpha=0.5,size=1,aes(x=bp_position,y=cM,color=sex)) +
 #     theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
      #labs(x="",y="Genetic distance (cM)",color="Sex") +
      theme(strip.background =element_rect(fill="white")) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      scale_color_manual(values=c( 
        'darkgreen',
        'purple' #Blue #X
      )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none") 

    marey_A3<-ggplot(subset(reshaped_map, reshaped_map$LG=="A3")) + 
      geom_line(aes(x=bp_position,y=cM,color=sex), lineend="round", linejoin="round") +
      theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
      #geom_point(alpha=0.5,size=1,aes(x=bp_position,y=cM,color=sex)) +
   #   theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
      #labs(x="",y="Genetic distance (cM)",color="Sex") +
      theme(strip.background =element_rect(fill="white")) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      scale_color_manual(values=c( 
        'darkgreen',
        'purple' #Blue #X
      )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none") 
    
    marey_A4<-ggplot(subset(reshaped_map, reshaped_map$LG=="A4")) + 
      geom_line(aes(x=bp_position,y=cM,color=sex), lineend="round", linejoin="round") +
  theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
  #      geom_point(alpha=0.5,size=1,aes(x=bp_position,y=cM,color=sex)) +
#      theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
      #labs(x="",y="Genetic distance (cM)",color="Sex") +
      theme(strip.background =element_rect(fill="white")) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      scale_color_manual(values=c( 
        'darkgreen',
        'purple' #Blue #X
      )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none") 
    
    marey_XY<-ggplot(subset(reshaped_map, reshaped_map$LG=="XY")) + 
      geom_line(aes(x=bp_position,y=cM,color=sex), lineend="round", linejoin="round") +
      theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
      #geom_point(alpha=0.5,size=1,aes(x=bp_position,y=cM,color=sex)) +
    #  theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
      #labs(x="",y="Genetic distance (cM)",color="Sex") +
      theme(strip.background =element_rect(fill="white")) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      scale_color_manual(values=c( 
        'darkgreen',
        'purple' #Blue #X
      )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none") 
    

#    max((subset(reshaped_map, reshaped_map$LG=="A1"))$bp_position)
 #   max((subset(reshaped_map, reshaped_map$LG=="A2"))$bp_position)
  #  max((subset(reshaped_map, reshaped_map$LG=="A3"))$bp_position)
   # max((subset(reshaped_map, reshaped_map$LG=="A4"))$bp_position)
    #max((subset(reshaped_map, reshaped_map$LG=="XY"))$bp_position)
    
    
    complete_marey<-grid.arrange(marey_A1,marey_A2, marey_A3, marey_A4, marey_XY, ncol=5, widths=c(3.88306085,2.78276620,1.71702641,1.35195936,2.38962673))
    
    
    
    
    
    
    
   # ggsave(filename = "Fig_1_A_Marey_9-13.pdf", plot = egg::set_panel_size(p=marey_plot_m_f, width=unit(10, "cm"), height=unit(3, "cm")))
    
  
    ################ 5.1.2 Fig 1B Sex differences in recombination ################
    #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
    
    complete_windowed_dataset$crossovers_diff<-(complete_windowed_dataset$crossovers_male_sum-complete_windowed_dataset$crossovers_female_sum)
    
    colnames(complete_windowed_dataset)
    
    complete_windowed_dataset <- complete_windowed_dataset %>%
      mutate(
        crossovers_binary = case_when(
          ((crossovers_male_sum > crossovers_female_sum) == TRUE) ~ 1,
          ((crossovers_male_sum < crossovers_female_sum) == TRUE) ~ 0))
    
    #View(reshaped_map)
    
   # plotting_windowed_map<-complete_windowed_dataset %>% 
 #     dplyr::select(
#        LG , window_end.x ,crossovers_male_sum , crossovers_female_sum , crossovers_sex_averaged_sum , crossovers_diff , male_distortion_mean , female_distortion_mean , n_genes_Yunchen , sum_Length_Yunchen ,N_TEs_all , sum_length_TEs_all , N_DNA_TEs_all , sum_length_DNA_TEs_all , N_RNA_TEs_all , sum_length_RNA_TEs_all)
    
      plotting_windowed_map<-complete_windowed_dataset %>% 
        dplyr::select(
          LG , window_end.x , crossovers_diff, crossovers_binary )

    reshaped_windowed_map<-plotting_windowed_map%>%pivot_longer(!c(LG,window_end.x), names_to = "chr", values_to = "crossovers_diff")
    
   # reshaped_map<-fullmap%>%pivot_longer(c(cm_male, cm_female), names_to = "sex", values_to = "cM")
    
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG4"] <- "XY"
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG1"] <- "A1"
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG3"] <- "A2"
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG5"] <- "A3"
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG2"] <- "A4"
    reshaped_windowed_map$LG = factor(reshaped_windowed_map$LG, levels=c('A1','A2','A3','A4','XY'))
    
    
    
    crossover_diff_plot <-
      ggplot(reshaped_windowed_map) + 
      geom_line(alpha=1,aes(x=window_end.x,y=crossovers_diff))+
      facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
      theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
      labs(x="",y="Difference in\ncrossover number") +
      theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
     # scale_color_manual(values=c( 
    #    'darkgreen',
    #    'purple' #Blue #X
    #  )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines")) 
    #Saved as pdf 5"x10"
    #View(reshaped_windowed_map)
    
    max((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A1"))$crossovers_diff, na.rm=T)
    min((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A1"))$crossovers_diff, na.rm=T)
    max((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A2"))$crossovers_diff, na.rm=T)
    min((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A2"))$crossovers_diff, na.rm=T)
    max((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A3"))$crossovers_diff, na.rm=T)
    min((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A3"))$crossovers_diff, na.rm=T)
    max((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A4"))$crossovers_diff, na.rm=T)
    min((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A4"))$crossovers_diff, na.rm=T)
    max((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="XY"))$crossovers_diff, na.rm=T)
    min((subset(reshaped_windowed_map, reshaped_windowed_map$LG=="XY"))$crossovers_diff, na.rm=T)
    
    crossover_diff_plot_A1 <-
      ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A1")) + 
      geom_line(alpha=1,aes(x=window_end.x,y=crossovers_diff))+
      #facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
      theme_bw()  + theme_bw(base_size = 18) + ylim(-35,20)+
      labs(x="",y="Difference in\ncrossover number") +
      theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      # scale_color_manual(values=c( 
      #    'darkgreen',
      #    'purple' #Blue #X
      #  )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines")) 

    
    crossover_diff_plot_A2 <-
      ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A2")) + 
      geom_line(alpha=1,aes(x=window_end.x,y=crossovers_diff))+
      #facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
      theme_bw()  + theme_bw(base_size = 18) + ylim(-35,20)+
      labs(x="",y="Difference in\ncrossover number") +
      theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      # scale_color_manual(values=c( 
      #    'darkgreen',
      #    'purple' #Blue #X
      #  )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines")) 
    
    
    crossover_diff_plot_A3 <-
      ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A3")) + 
      geom_line(alpha=1,aes(x=window_end.x,y=crossovers_diff))+
      #facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
      theme_bw()  + theme_bw(base_size = 18) + ylim(-35,20)+
      labs(x="",y="Difference in\ncrossover number") +
      theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      # scale_color_manual(values=c( 
      #    'darkgreen',
      #    'purple' #Blue #X
      #  )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines")) 
    
    crossover_diff_plot_A4 <-
      ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG=="A4")) + 
      geom_line(alpha=1,aes(x=window_end.x,y=crossovers_diff))+
      #facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
      theme_bw()  + theme_bw(base_size = 18) + ylim(-35,20)+
      labs(x="",y="Difference in\ncrossover number") +
      theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      # scale_color_manual(values=c( 
      #    'darkgreen',
      #    'purple' #Blue #X
      #  )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))   
    
    crossover_diff_plot_XY <-
      ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG=="XY")) + 
      geom_line(alpha=1,aes(x=window_end.x,y=crossovers_diff))+
      #facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
      theme_bw()  + theme_bw(base_size = 18) + ylim(-35,20)+
      labs(x="",y="Difference in\ncrossover number") +
      theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      # scale_color_manual(values=c( 
      #    'darkgreen',
      #    'purple' #Blue #X
      #  )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines")) 
    
    complete_crossover_diff<-grid.arrange(crossover_diff_plot_A1,crossover_diff_plot_A2, crossover_diff_plot_A3, crossover_diff_plot_A4, crossover_diff_plot_XY, ncol=5, widths=c(3.88306085,2.78276620,1.71702641,1.35195936,2.38962673))
    
    
    
    
    ################ 5.1.3 Fig 1C Plot distortion ################
    
    
    plotting_map<-fullmap %>% dplyr::select(LG, bp_position, male_distortion,female_distortion)
    
    
    reshaped_map<-plotting_map%>%pivot_longer(!c(LG,bp_position), names_to = "chr", values_to = "distortion")
    
    reshaped_map<-plotting_map%>%pivot_longer(c(male_distortion, female_distortion), names_to = "sex", values_to = "distortion")

    
    
    reshaped_map<-subset(reshaped_map,!is.na(reshaped_map$distortion))
    
    reshaped_map$LG[reshaped_map$LG == "LG4"] <- "XY"
    reshaped_map$LG[reshaped_map$LG == "LG1"] <- "A1"
    reshaped_map$LG[reshaped_map$LG == "LG3"] <- "A2"
    reshaped_map$LG[reshaped_map$LG == "LG5"] <- "A3"
    reshaped_map$LG[reshaped_map$LG == "LG2"] <- "A4"
    reshaped_map$LG = factor(reshaped_map$LG, levels=c('A1','A2','A3','A4','XY'))
    
    
    
    
    #View(reshaped_map)

      distortion_plot_m_f <-
      ggplot(reshaped_map) + 
      geom_line(aes(x=bp_position,y=distortion,color=sex)) + 
      #geom_line(aes(x=bp_position,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
      facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
      theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
      labs(x="",y="Segregation\ndistortion (LOD)",color="Sex") +
      theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      scale_color_manual(values=c( 
        'darkgreen',
        'purple' #Blue #X
      )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")+ 
    #theme(axis.title.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))
    geom_hline(yintercept=6.634897, linetype="dashed", color = "red")+
      geom_hline(yintercept=3.841459, linetype="dotted", color = "red")
    
      
  
      
      distortion_plot_m_f_A1 <-
        ggplot(subset(reshaped_map, reshaped_map$LG=="A1")) + 
        geom_line(aes(x=bp_position,y=distortion,color=sex)) + 
        #geom_line(aes(x=bp_position,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
       # facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,15)+
        labs(x="",y="Segregation\ndistortion (LOD)",color="Sex") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgreen',
          'purple' #Blue #X
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")+ 
        #theme(axis.title.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))
        geom_hline(yintercept=6.634897, linetype="dashed", color = "red")+
        geom_hline(yintercept=3.841459, linetype="dotted", color = "red")
      
      
      distortion_plot_m_f_A2 <-
        ggplot(subset(reshaped_map, reshaped_map$LG=="A2")) + 
        geom_line(aes(x=bp_position,y=distortion,color=sex)) + 
        #geom_line(aes(x=bp_position,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        # facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,15)+
        labs(x="",y="Segregation\ndistortion (LOD)",color="Sex") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgreen',
          'purple' #Blue #X
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")+ 
        #theme(axis.title.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))
        geom_hline(yintercept=6.634897, linetype="dashed", color = "red")+
        geom_hline(yintercept=3.841459, linetype="dotted", color = "red")
      
      distortion_plot_m_f_A3 <-
        ggplot(subset(reshaped_map, reshaped_map$LG=="A3")) + 
        geom_line(aes(x=bp_position,y=distortion,color=sex)) + 
        #geom_line(aes(x=bp_position,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        # facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,15)+
        labs(x="",y="Segregation\ndistortion (LOD)",color="Sex") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgreen',
          'purple' #Blue #X
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")+ 
        #theme(axis.title.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))
        geom_hline(yintercept=6.634897, linetype="dashed", color = "red")+
        geom_hline(yintercept=3.841459, linetype="dotted", color = "red")
      
      distortion_plot_m_f_A4 <-
        ggplot(subset(reshaped_map, reshaped_map$LG=="A4")) + 
        geom_line(aes(x=bp_position,y=distortion,color=sex)) + 
        #geom_line(aes(x=bp_position,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        # facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,15)+
        labs(x="",y="Segregation\ndistortion (LOD)",color="Sex") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgreen',
          'purple' #Blue #X
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")+ 
        #theme(axis.title.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))
        geom_hline(yintercept=6.634897, linetype="dashed", color = "red")+
        geom_hline(yintercept=3.841459, linetype="dotted", color = "red")
      
      distortion_plot_m_f_XY <-
        ggplot(subset(reshaped_map, reshaped_map$LG=="XY")) + 
        geom_line(aes(x=bp_position,y=distortion,color=sex)) + 
        #geom_line(aes(x=bp_position,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        # facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,15)+
        labs(x="",y="Segregation\ndistortion (LOD)",color="Sex") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgreen',
          'purple' #Blue #X
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")+ 
        #theme(axis.title.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))
        geom_hline(yintercept=6.634897, linetype="dashed", color = "red")+
        geom_hline(yintercept=3.841459, linetype="dotted", color = "red")
      
      
      
      complete_distortion_plot<-grid.arrange(distortion_plot_m_f_A1,distortion_plot_m_f_A2, distortion_plot_m_f_A3, distortion_plot_m_f_A4, distortion_plot_m_f_XY, ncol=5, widths=c(3.88306085,2.78276620,1.71702641,1.35195936,2.38962673))
      
      
      
      
      
      
    
    ################ 5.1.4 Fig 1D Plot gene density ################

    plotting_windowed_map<-complete_windowed_dataset %>% 
      dplyr::select(
        LG , window_end.x , n_genes_Yunchen )
    
    reshaped_windowed_map<-plotting_windowed_map%>%pivot_longer(!c(LG,window_end.x), names_to = "chr", values_to = "n_genes")
    
    # reshaped_map<-fullmap%>%pivot_longer(c(cm_male, cm_female), names_to = "sex", values_to = "cM")
    
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG4"] <- "XY"
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG1"] <- "A1"
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG3"] <- "A2"
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG5"] <- "A3"
    reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG2"] <- "A4"
    reshaped_windowed_map$LG = factor(reshaped_windowed_map$LG, levels=c('A1','A2','A3','A4','XY'))
    
    
   # View(reshaped_windowed_map)
    reshaped_map<-subset(reshaped_map,!is.na(reshaped_map$distortion))
    unique(reshaped_windowed_map$LG)
    
      gene_density_plot <-
      ggplot(reshaped_windowed_map) + 
      geom_line(aes(x=window_end.x,y=n_genes))+
      facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
      theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
      labs(x="",y="Genes per Mb") +
      theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
      # scale_color_manual(values=c( 
      #    'darkgreen',
      #    'purple' #Blue #X
      #  )) + 
      theme(axis.title.x=element_blank(), axis.text.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines")) 
    #Saved as pdf 5"x10"
    
      
      
      gene_density_plot_A1 <-
        ggplot(subset(reshaped_windowed_map,reshaped_windowed_map$LG=="A1")) + 
        geom_line(aes(x=window_end.x,y=n_genes))+
   #     facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
        labs(x="",y="Genes per Mb") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        # scale_color_manual(values=c( 
        #    'darkgreen',
        #    'purple' #Blue #X
        #  )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines")) 
      #Saved as pdf 5"x10"
      
      
      
      gene_density_plot_A2 <-
        ggplot(subset(reshaped_windowed_map,reshaped_windowed_map$LG=="A2")) + 
        geom_line(aes(x=window_end.x,y=n_genes))+
        #     facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
        labs(x="",y="Genes per Mb") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        # scale_color_manual(values=c( 
        #    'darkgreen',
        #    'purple' #Blue #X
        #  )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines")) 
      #Saved as pdf 5"x10"
      
      gene_density_plot_A3 <-
        ggplot(subset(reshaped_windowed_map,reshaped_windowed_map$LG=="A3")) + 
        geom_line(aes(x=window_end.x,y=n_genes))+
        #     facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
        labs(x="",y="Genes per Mb") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        # scale_color_manual(values=c( 
        #    'darkgreen',
        #    'purple' #Blue #X
        #  )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines")) 
      #Saved as pdf 5"x10"
      
      gene_density_plot_A4 <-
        ggplot(subset(reshaped_windowed_map,reshaped_windowed_map$LG=="A4")) + 
        geom_line(aes(x=window_end.x,y=n_genes))+
        #     facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
        labs(x="",y="Genes per Mb") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        # scale_color_manual(values=c( 
        #    'darkgreen',
        #    'purple' #Blue #X
        #  )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines")) 
      #Saved as pdf 5"x10"
      
      gene_density_plot_XY <-
        ggplot(subset(reshaped_windowed_map,reshaped_windowed_map$LG=="XY")) + 
        geom_line(aes(x=window_end.x,y=n_genes))+
        #     facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,120)+
        labs(x="",y="Genes per Mb") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        # scale_color_manual(values=c( 
        #    'darkgreen',
        #    'purple' #Blue #X
        #  )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines")) 
      #Saved as pdf 5"x10"
    
    
      complete_gene_density<-grid.arrange(gene_density_plot_A1,gene_density_plot_A2, gene_density_plot_A3, gene_density_plot_A4, gene_density_plot_XY, ncol=5, widths=c(3.88306085,2.78276620,1.71702641,1.35195936,2.38962673))
      
      
      
      
      
      
    ################ 5.1.5 Fig 1E TE density ################
    
      
      plotting_windowed_map<-complete_windowed_dataset %>% 
        dplyr::select(
          LG , window_end.x , N_DNA_TEs_all, N_RNA_TEs_all)
      
      reshaped_windowed_map<-plotting_windowed_map%>%pivot_longer(!c(LG,window_end.x), names_to = "TE_type", values_to = "n_TEs")
      
      reshaped_map<-plotting_windowed_map%>%pivot_longer(c(N_DNA_TEs_all, N_RNA_TEs_all), names_to = "TE_type", values_to = "count")
      
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG4"] <- "XY"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG1"] <- "A1"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG3"] <- "A2"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG5"] <- "A3"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG2"] <- "A4"
      reshaped_windowed_map$LG = factor(reshaped_windowed_map$LG, levels=c('A1','A2','A3','A4','XY'))
      
      
     # View(reshaped_windowed_map)
      reshaped_windowed_map$Mb<-reshaped_windowed_map$window_end.x/1000000
      
      TE_plot <-
        ggplot(reshaped_windowed_map) + 
        geom_line(alpha=0.5,aes(x=Mb,y=n_TEs,color=TE_type)) + 
        #  geom_line(aes(x=mb,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
        labs(x="",y="TEs per mb",color="TE_type") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgoldenrod2', #DNA TEs
          'cyan4' #RNA TEs
        )) + 
        theme(axis.title.x=element_blank())+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")
        

      
      TE_plot_A1 <-
        ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG =="A1")) + 
        geom_line(alpha=0.5,aes(x=Mb,y=n_TEs,color=TE_type)) + 
        #  geom_line(aes(x=mb,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
      #  facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,800)+
        labs(x="",y="TEs per mb",color="TE_type") +
        theme(strip.background =element_rect(fill="white"))  +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgoldenrod2', #DNA TEs
          'cyan4' #RNA TEs
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")
      
      
      TE_plot_A2 <-
        ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG =="A2")) + 
        geom_line(alpha=0.5,aes(x=Mb,y=n_TEs,color=TE_type)) + 
        #  geom_line(aes(x=mb,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        #  facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) +  ylim(0,800) +
        labs(x="",y="TEs per mb",color="TE_type") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgoldenrod2', #DNA TEs
          'cyan4' #RNA TEs
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")
      
      TE_plot_A3 <-
        ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG =="A3")) + 
        geom_line(alpha=0.5,aes(x=Mb,y=n_TEs,color=TE_type)) + 
        #  geom_line(aes(x=mb,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        #  facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,800) +
        labs(x="",y="TEs per mb",color="TE_type") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgoldenrod2', #DNA TEs
          'cyan4' #RNA TEs
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")
      
      TE_plot_A4 <-
        ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG =="A4")) + 
        geom_line(alpha=0.5,aes(x=Mb,y=n_TEs,color=TE_type)) + 
        #  geom_line(aes(x=mb,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        #  facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,800) +
        labs(x="",y="TEs per mb",color="TE_type") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgoldenrod2', #DNA TEs
          'cyan4' #RNA TEs
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")
      
      TE_plot_XY <-
        ggplot(subset(reshaped_windowed_map, reshaped_windowed_map$LG =="XY")) + 
        geom_line(alpha=0.5,aes(x=Mb,y=n_TEs,color=TE_type)) + 
        #  geom_line(aes(x=mb,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        #  facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + ylim(0,800) +
        labs(x="",y="TEs per mb",color="TE_type") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgoldenrod2', #DNA TEs
          'cyan4' #RNA TEs
        )) + 
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.title.y=element_blank(), axis.text.y=element_blank() )+ theme(panel.spacing = unit(0.1, "lines"))+ theme(legend.position = "none")
      
      
      
      
      complete_TE_density<-grid.arrange(TE_plot_A1,TE_plot_A2, TE_plot_A3, TE_plot_A4, TE_plot_XY, ncol=5, widths=c(3.88306085,2.78276620,1.71702641,1.35195936,2.38962673))
      
      
      
      
      
      ################ 5.1.6 Combine Fig 1 ################
      
#      grid.arrange(marey_plot_m_f, crossover_diff_plot, distortion_plot_m_f, gene_density_plot, TE_plot, nrow=5)
 
      
      y_marey <- gtable_filter(ggplotGrob(marey_A1), "guide-box") 
      # grid.draw(legend)    # Make sure the legend has been extracted
      
      y_marey <- gtable_filter(ggplotGrob(marey_A1), "ylab-l") 
      grid.draw(y_marey)
      ylab-l
           
      combined_figure_1<-grid.arrange(complete_marey, complete_crossover_diff, complete_distortion_plot, complete_gene_density,complete_TE_density, nrow=5)

      
      ################ 5.1.7 Supplementary chromosome figures ################
      
      ################ 5.1.7.1 S1 Fixed differences between X and Y in population data ################
      
      fixed_diffs<-read.table("Proportion_fixed_XY.txt")
      colnames(fixed_diffs)<-c("Mb", "Prop_fixed")
      
      fixed_diffs_plot <-
        ggplot(fixed_diffs) + 
        geom_line(alpha=1,aes(x=Mb,y=Prop_fixed)) + 
        #  geom_line(aes(x=mb,y=CM_p),linetype="dashed",span=5,se=FALSE,color="#ffb14e")+
        #  facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18)  +
        labs(x="Position on sex chromosome (Mb)",y="Proportion fixed\nX-Y differences") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
      
      #ggsave("Figures/Proportion_fixed_XY.pdf", width = 8, height = 4, units = "in")
      #ggsave("Figures/Proportion_fixed_XY.png", width = 800, height = 400, units = "px", dpi=72)    
      
      
      ################ 5.1.7.6 S2 Genotype frequencies ################
      #Code is also here to just display raw genotypes, only needs a little rearranging
      
      genotypes<-read.csv("Joanna_final/phys/all_geno.csv")
      head(genotypes)
      genotypes<-genotypes  %>% separate(., CHR, into = c("LG", "Position"), sep = "\\_") 
      genotypes<-genotypes %>% mutate(n_11 = rowSums(.[3:188]=="11"))
      genotypes<-genotypes %>% mutate(n_12 = rowSums(.[3:188]=="12"))
      genotypes<-genotypes %>% mutate(n_21 = rowSums(.[3:188]=="21"))
      genotypes<-genotypes %>% mutate(n_22 = rowSums(.[3:188]=="22"))
      head(genotypes[189:192])
      
      chromosome_gts<-subset(genotypes, genotypes$LG=="LG5")
      plot(chromosome_gts$Position, chromosome_gts$n_11)
      plot(chromosome_gts$Position, chromosome_gts$n_12)
      plot(chromosome_gts$Position, chromosome_gts$n_22)
      
      #11 = paternal allele, 22 = maternal allele
      
      plotting_genotypes<-dplyr::select(genotypes, c(LG, Position, n_11,n_12,n_21,n_22))
      
      # plotting_genotypes<-read.csv("Genotypes_for_all_markers.csv")
      plotting_genotypes$prop_11<-(plotting_genotypes$n_11/(plotting_genotypes$n_11+plotting_genotypes$n_12+plotting_genotypes$n_21+plotting_genotypes$n_22))
      
      plotting_genotypes<-plotting_genotypes%>%rowwise()%>%mutate(., prop_11=(n_11/sum(n_11,n_12,n_21,n_22)))
      plotting_genotypes<-plotting_genotypes%>%rowwise()%>%mutate(., prop_het=((n_12+n_21)/sum(n_11,n_12,n_21,n_22)))
      plotting_genotypes<-plotting_genotypes%>%rowwise()%>%mutate(., prop_22=(n_22/sum(n_11,n_12,n_21,n_22)))
      
      plotting_genotypes<-plotting_genotypes%>%rowwise()%>%mutate(., prop_1_dad=(
        (n_11+n_12)/
          (sum(n_11,n_12,n_21,n_22)))
        )
      plotting_genotypes<-plotting_genotypes%>%rowwise()%>%mutate(., prop_2_dad=(
        (n_22+n_21)/
        (sum(n_11,n_12,n_21,n_22)))
      )
      plotting_genotypes<-plotting_genotypes%>%rowwise()%>%mutate(., prop_1_mum=(
        (n_11+n_21)/
          (sum(n_11,n_12,n_21,n_22)))
        )
      plotting_genotypes<-plotting_genotypes%>%rowwise()%>%mutate(., prop_2_mum=(
        (n_22+n_12)/
          (sum(n_11,n_12,n_21,n_22)))
        )
      
      
      
      View(plotting_genotypes)
      plotting_genotypes<-dplyr::select(plotting_genotypes, c(LG, Position, prop_11,prop_het,prop_22))
      
      plotting_genotypes<-dplyr::select(plotting_genotypes, c(LG, Position, prop_1_dad,prop_2_dad,prop_1_mum,prop_2_mum))
      
      
      reshaped_plotting_genotypes<-plotting_genotypes%>%pivot_longer(!c(LG,Position), names_to = "what_GT", values_to = "GT")

      
            
      reshaped_plotting_genotypes$LG[reshaped_plotting_genotypes$LG == "LG4"] <- "XY"
      reshaped_plotting_genotypes$LG[reshaped_plotting_genotypes$LG == "LG1"] <- "A1"
      reshaped_plotting_genotypes$LG[reshaped_plotting_genotypes$LG == "LG3"] <- "A2"
      reshaped_plotting_genotypes$LG[reshaped_plotting_genotypes$LG == "LG5"] <- "A3"
      reshaped_plotting_genotypes$LG[reshaped_plotting_genotypes$LG == "LG2"] <- "A4"
      reshaped_plotting_genotypes$LG = factor(reshaped_plotting_genotypes$LG, levels=c('A1','A2','A3','A4','XY'))
      reshaped_plotting_genotypes$Mb<-reshaped_plotting_genotypes$Position/1000000
      
      
      gt_plot <-
        ggplot(reshaped_plotting_genotypes) + 
        geom_line(aes(x=Mb,y=GT, color=what_GT, linetype=what_GT))+
        facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
        labs(x="",y="Proportion alleles passed on",color="Genotype") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +  #ylim(0,1) +
        scale_linetype_manual(values = c(rep("solid", 2), rep("dashed", 2))) +
        scale_color_manual(labels = c("AA", "BB","AB"),values=c( 
          'purple',
          'darkgreen',
          'purple',
          'darkgreen'
          #     'cyan4',
      #    '#C01B3A'
        )) + 
      #  # theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
        theme(panel.spacing = unit(0.1, "lines")) + theme(legend.position = "none")

            
      # write.csv(plotting_genotypes, "Genotypes_for_all_markers.csv", quote=F, row.name=F)
      
      #ggsave("Figures/genotypes.pdf", width = 10, height = 4, units = "in")
      #ggsave("Figures/genotypes.png", width = 1000, height = 400, units = "px", dpi=72)    
      
      
      ################ 5.1.7.2 S2 Leaf and flower expression normalized to number of genes in window ################
      
      #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
      
      
      complete_windowed_dataset$expressed_leaf_prop<-complete_windowed_dataset$expressed_leaf_sum/complete_windowed_dataset$n_genes_Yunchen
      
      complete_windowed_dataset$expressed_flower_prop<-complete_windowed_dataset$expressed_flower_sum/complete_windowed_dataset$n_genes_Yunchen
      
      plotting_windowed_map<-complete_windowed_dataset %>% 
        dplyr::select(
          LG , window_end.x , expressed_leaf_prop, expressed_flower_prop )
      
      reshaped_windowed_map<-plotting_windowed_map%>%pivot_longer(!c(LG,window_end.x), names_to = "tissue", values_to = "expression")
      
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG4"] <- "XY"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG1"] <- "A1"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG3"] <- "A2"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG5"] <- "A3"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG2"] <- "A4"
      reshaped_windowed_map$LG = factor(reshaped_windowed_map$LG, levels=c('A1','A2','A3','A4','XY'))
      
      reshaped_windowed_map<-subset(reshaped_windowed_map, !is.na(reshaped_windowed_map$LG))
      reshaped_windowed_map$Mb<-reshaped_windowed_map$window_end.x/1000000
      
      normalized_expression_plot <-
        ggplot(reshaped_windowed_map) + 
        geom_line(alpha=.5,aes(x=Mb,y=expression, color=tissue))+
        facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
        labs(x="",y="Proportion leaf- or flower-\nexpressed per window",color="Tissue") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        # scale_color_manual(values=c( 
        #    'darkgreen',
        #    'purple' #Blue #X
        #  )) + 
       # theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
        theme(panel.spacing = unit(0.1, "lines")) #+ theme(legend.position = "none")
      #Cyan = proportion expressed in flower
      #coral = proportan expressed in leaf
      
      ################ 5.1.7.3 S3 Male- and female-biased leaf expression normalized to number of genes in window ################
      
      
      #complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
      
      
      complete_windowed_dataset$female_leaf_prop<-complete_windowed_dataset$leaf_female_bias_sum /complete_windowed_dataset$n_genes_Yunchen
      
      complete_windowed_dataset$male_leaf_prop<-complete_windowed_dataset$leaf_male_bias_sum/complete_windowed_dataset$n_genes_Yunchen
      
      plotting_windowed_map<-complete_windowed_dataset %>% 
        dplyr::select(
          LG , window_end.x , female_leaf_prop, male_leaf_prop )
      
      reshaped_windowed_map<-plotting_windowed_map%>%pivot_longer(!c(LG,window_end.x), names_to = "sex_bias", values_to = "expression")
      
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG4"] <- "XY"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG1"] <- "A1"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG3"] <- "A2"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG5"] <- "A3"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG2"] <- "A4"
      reshaped_windowed_map$LG = factor(reshaped_windowed_map$LG, levels=c('A1','A2','A3','A4','XY'))
      
      reshaped_windowed_map<-subset(reshaped_windowed_map, !is.na(reshaped_windowed_map$LG))
      reshaped_windowed_map$Mb<-reshaped_windowed_map$window_end.x/1000000
      
      normalized_sex_bias_leaf_plot <-
        ggplot(reshaped_windowed_map) + 
        geom_line(alpha=.5,aes(x=Mb,y=expression, color=sex_bias))+
        facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
        labs(x="",y="Proportion of genes with\nmale- or female-biased leaf expression",color="Sex") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
         scale_color_manual(values=c( 
            'darkgreen',
            'purple' #Blue #X
          )) + 
        # theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
        theme(panel.spacing = unit(0.1, "lines")) + theme(legend.position = "none")
      #purple = male bias
      #green = female bias
      
      
      ################ 5.1.7.4 S4 Male- and female-biased flower expression normalized to number of genes in window ################
      
      
      #complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
      
      
      complete_windowed_dataset$female_flower_prop<-complete_windowed_dataset$flower_female_bias_sum /complete_windowed_dataset$n_genes_Yunchen
      
      complete_windowed_dataset$male_flower_prop<-complete_windowed_dataset$flower_male_bias_sum/complete_windowed_dataset$n_genes_Yunchen
      
      plotting_windowed_map<-complete_windowed_dataset %>% 
        dplyr::select(
          LG , window_end.x , female_flower_prop, male_flower_prop )
      
      reshaped_windowed_map<-plotting_windowed_map%>%pivot_longer(!c(LG,window_end.x), names_to = "sex_bias", values_to = "expression")
      
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG4"] <- "XY"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG1"] <- "A1"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG3"] <- "A2"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG5"] <- "A3"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG2"] <- "A4"
      reshaped_windowed_map$LG = factor(reshaped_windowed_map$LG, levels=c('A1','A2','A3','A4','XY'))
      
      reshaped_windowed_map<-subset(reshaped_windowed_map, !is.na(reshaped_windowed_map$LG))
      reshaped_windowed_map$Mb<-reshaped_windowed_map$window_end.x/1000000
      
      normalized_sex_bias_flower_plot <-
        ggplot(reshaped_windowed_map) + 
        geom_line(alpha=.5,aes(x=Mb,y=expression, color=sex_bias))+
        facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
        labs(x="",y="Proportion of genes with\nmale- or female-biased flower expression",color="Sex") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgreen',
          'purple' #Blue #X
        )) + 
        # theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
        theme(panel.spacing = unit(0.1, "lines")) + theme(legend.position = "none")
      #purple = male bias
      #green = female bias
      
      
      
      
      
      ################ 5.1.7.5 S5 Pollen- and pollen-tube-biased expression normalized to number of genes in window ################
      
      
      #complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
      
      complete_windowed_dataset$pollen_biased_sum_prop<-complete_windowed_dataset$pollen_biased_sum/complete_windowed_dataset$n_genes_Yunchen
      
      complete_windowed_dataset$pollen_tube_biased_prop<-complete_windowed_dataset$pollen_tube_biased_sum/complete_windowed_dataset$n_genes_Yunchen
      
      plotting_windowed_map<-complete_windowed_dataset %>% 
        dplyr::select(
          LG , window_end.x , pollen_biased_sum_prop, pollen_tube_biased_prop )
      
      reshaped_windowed_map<-plotting_windowed_map%>%pivot_longer(!c(LG,window_end.x), names_to = "tissue", values_to = "expression")
      
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG4"] <- "XY"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG1"] <- "A1"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG3"] <- "A2"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG5"] <- "A3"
      reshaped_windowed_map$LG[reshaped_windowed_map$LG == "LG2"] <- "A4"
      reshaped_windowed_map$LG = factor(reshaped_windowed_map$LG, levels=c('A1','A2','A3','A4','XY'))
      
      reshaped_windowed_map<-subset(reshaped_windowed_map, !is.na(reshaped_windowed_map$LG))
      reshaped_windowed_map$Mb<-reshaped_windowed_map$window_end.x/1000000
      
      normalized_pollen_bias_plot <-
        ggplot(reshaped_windowed_map) + 
        geom_line(alpha=.5,aes(x=Mb,y=expression, color=tissue))+
        facet_grid(. ~ LG, scales = 'free_x', space = 'free_x') +
        theme_bw()  + theme_bw(base_size = 18) + #xlim(0,20)+
        labs(x="",y="Proportion of genes with pollen- or\npollen tube-biased flower expression",color="Sex") +
        theme(strip.background =element_rect(fill="white")) + #xlim(120,125) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
        scale_color_manual(values=c( 
          'darkgoldenrod2',
          'cyan4'
        )) + 
        # theme(axis.title.x=element_blank(), axis.text.x=element_blank())+
        theme(panel.spacing = unit(0.1, "lines")) + theme(legend.position = "none")
      #purple = male bias
      #green = female bias
      
      
      ################ 5.2 Figure 2 ################
      
      ################ 5.2.1 Add functions ################
      
      
            #Plots based on https://towardsdatascience.com/customizable-correlation-plots-in-r-b1d2856a4b05
      
      cors <- function(df) { 
        # turn all three matrices (r, n, and P into a data frame)
        M <- Hmisc::rcorr(as.matrix(df))
        # return the three data frames in a list return(Mdf)
        Mdf <- map(M, ~data.frame(.x))
      }
      
      formatted_cors <- function(df){
        cors(df) %>%
          map(~rownames_to_column(.x, var="measure1")) %>%
          map(~pivot_longer(.x, -measure1, "measure2")) %>% 
          bind_rows(.id = "id") %>%
          pivot_wider(names_from = id, values_from = value) %>%
          mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), r_if_sig = ifelse(P <.05, r, NA)) 
      }
      
   #   View(formatted_cors(mtcars) )
      
 
      ################ 5.2.2 Whole-genome correlation plot ################

      #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
      complete_windowed_dataset$crossovers_diff<-(complete_windowed_dataset$crossovers_male_sum-complete_windowed_dataset$crossovers_female_sum)
      
      
      correlation_subset<-complete_windowed_dataset%>%dplyr::select(LG, n_genes_Yunchen , N_DNA_TEs_all , N_RNA_TEs_all , male_distortion_mean , female_distortion_mean , Mendelian_distortion_mean , crossovers_male_sum , crossovers_female_sum , crossovers_sex_averaged_sum , crossovers_diff , expressed_leaf_sum ,  leaf_male_bias_sum , leaf_female_bias_sum ,expressed_flower_sum , flower_male_bias_sum , flower_female_bias_sum , pollen_biased_sum , pollen_tube_biased_sum )

      whole_genome_corr_subset<-subset(correlation_subset, correlation_subset$LG %in% c("LG1", "LG2", "LG3", "LG4", "LG5"))
    #  whole_genome_corr_subset<-subset(correlation_subset, correlation_subset$LG %in% c("A1", "A2", "A3", "A4", "XY"))
      whole_genome_corr_subset<-whole_genome_corr_subset %>% dplyr::select(!LG)
      
      colnames(whole_genome_corr_subset)
     
      # colnames(whole_genome_corr_subset)<-c(
       # "Genes", "DNA_TEs", "RNA_TEs", 
      #  "M_TRD", "F_TRD", "TRD", 
       # "XO_M", "XO_F", "XO", "XO_diff", 
      #  "Leaf_exp", "Leaf_M", "Leaf_F", 
      #  "Flower_exp", "Flower_M", "Flower_F", 
      #  "Pollen", "Pol_tube")
      
      colnames(whole_genome_corr_subset)<-c(
        "Genes per window", "DNA TEs per window", "RNA TEs per window", 
        "Male distortion", "Female distortion", "Mendelian distortion", 
        "Male crossovers", "Female crossovers", "Total crossovers", "Crossover difference", 
        "Leaf expressed genes", "Leaf male biased", "Leaf female biased", 
        "Flower expressed genes", "Flower male biased", "Flower female biased", 
        "Pollen biased", "Pollen tube biased")
      
)
      formatted_cors(whole_genome_corr_subset) %>%
        mutate(measure2a=str_replace_all(measure2, "\\.", " ")) %>%
        mutate(measure1 = factor(measure1, levels=c(
       #   "Genes per window",  "DNA_TEs",  "RNA_TEs",  "M_TRD",  "F_TRD", "TRD",  "Leaf_exp",  "Leaf_M", "Leaf_F",  "Flower_exp",  "Flower_M",  "Flower_F",  "Pollen","Pol_tube"
          "Genes per window", "DNA TEs per window", "RNA TEs per window", 
          "Male distortion", "Female distortion", "Mendelian distortion", "Leaf expressed genes", "Leaf male biased", "Leaf female biased", 
          "Flower expressed genes", "Flower male biased", "Flower female biased", 
          "Pollen biased", "Pollen tube biased"
          ))) %>%
        mutate(measure2 = factor(measure2a, levels=c(
          "Genes per window",  "Male crossovers", "Female crossovers", "Total crossovers", "Crossover difference" #"Genes per window", "XO_diff",  "XO",  "XO_F",  "XO_M"
          ))) %>%  drop_na(., c(measure1,measure2))  %>%   ## to get the rect filled %>% 
        ggplot(aes(measure1, measure2, col=r)) +
        geom_tile(col="black", fill="white") + geom_point(aes(size = abs(r)), shape=15) + labs(x = NULL, y = NULL, col = "Pearson's\nCorrelation", title=paste("Genome-wide windowed correlations")) + theme_classic() +
        scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        #          scale_size(range=c(1,11), guide=NULL) + 
        scale_size(range=c(1,11), guide=NULL) + 
        theme(axis.text.x = element_text(angle=90))  

      
      correlation_output<-formatted_cors(whole_genome_corr_subset) %>%
        mutate(measure2a=str_replace_all(measure2, "\\.", " "))
      
      correlation_output$measure2<-correlation_output$measure2a
      correlation_output<-dplyr::select(correlation_output, -measure2a)
      
      write.csv(correlation_output, "Genome_wide_windowed_correlations.csv", quote=F, row.names=F)
      
      
      View(correlation_output)
      
      
      
      
    View ( formatted_cors(whole_genome_corr_subset) %>%
        mutate(measure1 = factor(measure1, levels=c(
          "Pol_tube",  "Pollen",  "Flower_F",  "Flower_M",  "Flower_exp", "Leaf_F",  "Leaf_M",  "Leaf_exp", "TRD",  "F_TRD",  "M_TRD",  "RNA_TEs",  "DNA_TEs", "Genes"
        ))) %>%
        mutate(measure2 = factor(measure2, levels=c(
          "XO_diff",  "XO",  "XO_F",  "XO_M"
        ))) %>% drop_na(., c(measure1,measure2)))
      
      
      
      formatted_cors(whole_genome_corr_subset) %>%
        mutate(measure1 = factor(measure1, levels=c("Pol_tube",  "Pollen",  "Flower_F",  "Flower_M",  "Flower_exp", "Leaf_F",  "Leaf_M",  "Leaf_exp",  "XO_diff",  "XO",  "XO_F",  "XO_M",  "TRD",  "F_TRD",  "M_TRD",  "RNA_TEs",  "DNA_TEs", "Genes"))) %>%
        mutate(measure2 = factor(measure2, levels=c("Pol_tube",  "Pollen",  "Flower_F",  "Flower_M",  "Flower_exp", "Leaf_F",  "Leaf_M",  "Leaf_exp",  "XO_diff",  "XO",  "XO_F",  "XO_M",  "TRD",  "F_TRD",  "M_TRD",  "RNA_TEs",  "DNA_TEs", "Genes"))) %>%       ## to get the rect filled 
        ggplot(aes(measure1, measure2, col=r)) +
        geom_tile(col="black", fill="white") + geom_point(aes(size = abs(r)), shape=15) + labs(x = NULL, y = NULL, col = "Pearson's\nCorrelation", title=paste("Genome-wide windowed correlations")) + theme_classic() +
        scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        #          scale_size(range=c(1,11), guide=NULL) + 
        scale_size(range=c(1,11), guide=NULL) + 
        theme(axis.text.x = element_text(angle=90))  
      
      formatted_cors_for_saving<-formatted_cors(whole_genome_corr_subset)        
      View(formatted_cors_for_saving)
      
      write.csv(formatted_cors_for_saving, "E:/Users/Joanna/Dropbox/Professional/University_of_Toronto/Genomics/LinkageMap/Replacement_sex_specific_recombination_map/Figure_correlations_9-13.csv")
      
      
      
      
      
      
      
      ################ 5.2.3 Single-chromosome correlation plot ################
      
      #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
      complete_windowed_dataset$crossovers_diff<-(complete_windowed_dataset$crossovers_male_sum-complete_windowed_dataset$crossovers_female_sum)
      
      
      correlation_subset<-complete_windowed_dataset%>%dplyr::select(LG, n_genes_Yunchen , N_DNA_TEs_all , N_RNA_TEs_all , male_distortion_mean , female_distortion_mean , Mendelian_distortion_mean , crossovers_male_sum , crossovers_female_sum , crossovers_sex_averaged_sum , crossovers_diff , expressed_leaf_sum ,  leaf_male_bias_sum , leaf_female_bias_sum ,expressed_flower_sum , flower_male_bias_sum , flower_female_bias_sum , pollen_biased_sum , pollen_tube_biased_sum )
      
      correlation_subset$LG[correlation_subset$LG == "LG4"] <- "XY"
      correlation_subset$LG[correlation_subset$LG == "LG1"] <- "A1"
      correlation_subset$LG[correlation_subset$LG == "LG3"] <- "A2"
      correlation_subset$LG[correlation_subset$LG == "LG5"] <- "A3"
      correlation_subset$LG[correlation_subset$LG == "LG2"] <- "A4"
      #  correlation_subset$LG = factor(reshaped_map$LG, levels=c('A1','A2','A3','A4','XY'))
      
      
      chromosome_correlation_subset<-subset(correlation_subset, correlation_subset$LG=="A1")
      main_title<-"A. A1 windowed correlations"
      
      chromosome_correlation_subset<-subset(correlation_subset, correlation_subset$LG=="A2")
      main_title<-"B. A2 windowed correlations"
      
      chromosome_correlation_subset<-subset(correlation_subset, correlation_subset$LG=="A3")
      main_title<-"C. A3 windowed correlations"
      
      chromosome_correlation_subset<-subset(correlation_subset, correlation_subset$LG=="A4")
      main_title<-"D. A4 windowed correlations"
      
      chromosome_correlation_subset<-subset(correlation_subset, correlation_subset$LG=="XY")
      main_title<-"E. XY windowed correlations"

      
      
      chromosome_correlation_subset<-chromosome_correlation_subset %>% dplyr::select(!LG)
      
      colnames(chromosome_correlation_subset)<-c(
        "Genes per window", "DNA TEs per window", "RNA TEs per window", 
        "Male distortion", "Female distortion", "Mendelian distortion", 
        "Male crossovers", "Female crossovers", "Total crossovers", "Crossover difference", 
        "Leaf expressed genes", "Leaf male biased", "Leaf female biased", 
        "Flower expressed genes", "Flower male biased", "Flower female biased", 
        "Pollen biased", "Pollen tube biased")
      
      
      
      
      corr_plot<-formatted_cors(chromosome_correlation_subset) %>%
        mutate(measure2a=str_replace_all(measure2, "\\.", " ")) %>%
        mutate(measure1 = factor(measure1, levels=c(
          #   "Genes per window",  "DNA_TEs",  "RNA_TEs",  "M_TRD",  "F_TRD", "TRD",  "Leaf_exp",  "Leaf_M", "Leaf_F",  "Flower_exp",  "Flower_M",  "Flower_F",  "Pollen","Pol_tube"
          "Genes per window", "DNA TEs per window", "RNA TEs per window", 
          "Male distortion", "Female distortion", "Mendelian distortion", "Leaf expressed genes", "Leaf male biased", "Leaf female biased", 
          "Flower expressed genes", "Flower male biased", "Flower female biased", 
          "Pollen biased", "Pollen tube biased"
        ))) %>%
        mutate(measure2 = factor(measure2a, levels=c(
          "Genes per window",  "Male crossovers", "Female crossovers", "Total crossovers", "Crossover difference" #"Genes per window", "XO_diff",  "XO",  "XO_F",  "XO_M"
        ))) %>%  drop_na(., c(measure1,measure2))  %>%   ## to get the rect filled %>% 
        ggplot(aes(measure1, measure2, col=r)) +
        geom_tile(col="black", fill="white") + geom_point(aes(size = abs(r)), shape=15) + labs(x = NULL, y = NULL, col = "Pearson's\nCorrelation", title=paste(main_title)) + theme_classic() +
        scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
        scale_x_discrete(expand=c(0,0)) +
        scale_y_discrete(expand=c(0,0)) +
        #          scale_size(range=c(1,11), guide=NULL) + 
        scale_size(range=c(1,11), guide=NULL) + 
        theme(axis.text.x = element_text(angle=90))  
      
      #ggsave("Figures/A1_correlations.pdf", width = 8, height = 4, units = "in")
      #ggsave("Figures/A1_correlations.png", width = 800, height = 400, units = "px", dpi=72)
      
      #ggsave("Figures/A2_correlations.pdf", width = 8, height = 4, units = "in")
      #ggsave("Figures/A2_correlations.png", width = 800, height = 400, units = "px", dpi=72)
      
      #ggsave("Figures/A3_correlations.pdf", width = 8, height = 4, units = "in")
      #ggsave("Figures/A3_correlations.png", width = 800, height = 400, units = "px", dpi=72)
      
      #ggsave("Figures/A4_correlations.pdf", width = 8, height = 4, units = "in")
      #ggsave("Figures/A4_correlations.png", width = 800, height = 400, units = "px", dpi=72)
      
      #ggsave("Figures/XY_correlations.pdf", width = 8, height = 4, units = "in")
      #ggsave("Figures/XY_correlations.png", width = 800, height = 400, units = "px", dpi=72)
  
      correlation_output<-formatted_cors(chromosome_correlation_subset) %>%
        mutate(measure2a=str_replace_all(measure2, "\\.", " "))
      
      correlation_output$measure2<-correlation_output$measure2a
      correlation_output<-dplyr::select(correlation_output, -measure2a)
      
      #write.csv(correlation_output, "A1_windowed_correlations.csv", quote=F, row.names=F)
      #write.csv(correlation_output, "A2_windowed_correlations.csv", quote=F, row.names=F)
      #write.csv(correlation_output, "A3_windowed_correlations.csv", quote=F, row.names=F)
      #write.csv(correlation_output, "A4_windowed_correlations.csv", quote=F, row.names=F)
      #write.csv(correlation_output, "XY_windowed_correlations.csv", quote=F, row.names=F)
      
      
      
      
      ################ 5.2.4 Whole-genome partial correlation plot ################
      
      
      #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
      complete_windowed_dataset$crossovers_diff<-(complete_windowed_dataset$crossovers_male_sum-complete_windowed_dataset$crossovers_female_sum)
      
      
      correlation_subset<-complete_windowed_dataset%>%dplyr::select(LG, n_genes_Yunchen , N_DNA_TEs_all , N_RNA_TEs_all , male_distortion_mean , female_distortion_mean , Mendelian_distortion_mean , crossovers_male_sum , crossovers_female_sum , crossovers_sex_averaged_sum , crossovers_diff , expressed_leaf_sum ,  leaf_male_bias_sum , leaf_female_bias_sum ,expressed_flower_sum , flower_male_bias_sum , flower_female_bias_sum , pollen_biased_sum , pollen_tube_biased_sum )
      
      whole_genome_corr_subset<-subset(correlation_subset, correlation_subset$LG %in% c("LG1", "LG2", "LG3", "LG4", "LG5"))
      #  whole_genome_corr_subset<-subset(correlation_subset, correlation_subset$LG %in% c("A1", "A2", "A3", "A4", "XY"))
      whole_genome_corr_subset<-whole_genome_corr_subset %>% dplyr::select(!LG)
      
  
        partial_correlations_results<-data.frame(matrix(ncol = 7, nrow = 0))
        colnames(partial_correlations_results)<-c("control","var_1","var_2","estimate","p","statistic","n")
        
              
      for (i in 2:length(colnames(whole_genome_corr_subset))){
  #      print(i)
       # var_1<-i
        for (j in 2:length(colnames(whole_genome_corr_subset))){
          if (i!=j){ 
         # var_2<-j
          #print(var_1)
          print(colnames(whole_genome_corr_subset)[i])
          #print(var_2)
          print(colnames(whole_genome_corr_subset)[j])
         # pcor.test(whole_genome_corr_subset[,i],whole_genome_corr_subset[,j],whole_genome_corr_subset[,"n_genes_Yunchen"])
          df<-data.frame(whole_genome_corr_subset[,i],whole_genome_corr_subset[,j],whole_genome_corr_subset[,"n_genes_Yunchen"])
          df<-df%>%drop_na(.)
          colnames(df)<-c(colnames(whole_genome_corr_subset)[i],colnames(whole_genome_corr_subset)[j],"n_genes_Yunchen")
         print( head(df))
         pcor_test_output<-pcor.test(df[,1], df[,2], df[3])
         print(pcor_test_output)
         outrow<-c(
           "n_genes_Yunchen",
           (colnames(whole_genome_corr_subset)[i]),
           (colnames(whole_genome_corr_subset)[j]),
           pcor_test_output$estimate,
           pcor_test_output$p.value,
           pcor_test_output$statistic,
           pcor_test_output$n
)  
         outrow<-data.frame(t(unlist(outrow)))
         colnames(outrow)<-c("control","var_1","var_2","estimate","p","statistic","n")
         print(outrow)
         partial_correlations_results<-rbind(partial_correlations_results,outrow)
          }
         }
      }
   
        
     
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^N_DNA_TEs_all", replacement = "DNA TEs per window")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^N_RNA_TEs_all", replacement = "RNA TEs per window")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^male_distortion_mean", replacement = "Male distortion")         
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^female_distortion_mean", replacement = "Female distortion")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^Mendelian_distortion_mean", replacement = "Mendelian distortion")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^crossovers_male_sum", replacement = "Male crossovers")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^crossovers_female_sum", replacement = "Female crossovers")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^crossovers_sex_averaged_sum", replacement = "Total crossovers")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^crossovers_diff", replacement = "Crossover difference")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^expressed_leaf_sum", replacement = "Leaf expressed genes")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^leaf_male_bias_sum", replacement = "Leaf male biased")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^leaf_female_bias_sum", replacement = "Leaf female biased")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^expressed_flower_sum", replacement = "Flower expressed genes")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^flower_male_bias_sum", replacement = "Flower male biased")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^flower_female_bias_sum", replacement = "Flower female biased")
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^pollen_biased_sum", replacement = "Pollen biased") 
        partial_correlations_results[2:3] <- lapply(partial_correlations_results[2:3], gsub, pattern = "^pollen_tube_biased_sum", replacement = "Pollen tube biased")
          

        partial_correlations_results$estimate<-as.numeric(partial_correlations_results$estimate)
        
        partial_correlations_results %>%
          mutate(var_1 = factor(var_1, levels=c(
            "DNA TEs per window", "RNA TEs per window", 
            "Male distortion", "Female distortion", "Mendelian distortion", "Leaf expressed genes", "Leaf male biased", "Leaf female biased", 
            "Flower expressed genes", "Flower male biased", "Flower female biased", 
            "Pollen biased", "Pollen tube biased"
          ))) %>%
          mutate(var_2 = factor(var_2, levels=c(
            "Male crossovers", "Female crossovers", "Total crossovers", "Crossover difference" #"Genes per window", "XO_diff",  "XO",  "XO_F",  "XO_M"
          ))) %>%  drop_na(., c(var_1,var_2))  %>%   ## to get the rect filled %>% 
          ggplot(aes(var_1,var_2, col=estimate)) +
          geom_tile(col="black", fill="white") + geom_point(aes(size = abs(estimate)), shape=15) + labs(x = NULL, y = NULL, col = "Pearson's\nCorrelation", title=paste("Genome-wide windowed partial correlations")) + theme_classic() +
          scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
          scale_x_discrete(expand=c(0,0)) +
          scale_y_discrete(expand=c(0,0)) +
          #          scale_size(range=c(1,11), guide=NULL) + 
          scale_size(range=c(1,11), guide=NULL) + 
          theme(axis.text.x = element_text(angle=90))  
        
        
        #ggsave("Figures/Genome_wide_partial_correlations.pdf", width = 8, height = 4, units = "in")
        #ggsave("Figures/Genome_wide_partial_correlations.png", width = 800, height = 400, units = "px", dpi=72)    
        
        write.csv(partial_correlations_results, "Genome_wide_windowed_partial_correlations.csv", quote=F, row.names=F)
        
        
      
      
      ################ 5.2.5 Single-chromosome partial correlation plots ################
        
       
        #    complete_windowed_dataset<-read.csv("Windowed_analyses/complete_windowed_data_3.3_8-30-2021.csv")
        complete_windowed_dataset$crossovers_diff<-(complete_windowed_dataset$crossovers_male_sum-complete_windowed_dataset$crossovers_female_sum)
        
        
        whole_genome_corr_subset<-complete_windowed_dataset%>%dplyr::select(LG, n_genes_Yunchen , N_DNA_TEs_all , N_RNA_TEs_all , male_distortion_mean , female_distortion_mean , Mendelian_distortion_mean , crossovers_male_sum , crossovers_female_sum , crossovers_sex_averaged_sum , crossovers_diff , expressed_leaf_sum ,  leaf_male_bias_sum , leaf_female_bias_sum ,expressed_flower_sum , flower_male_bias_sum , flower_female_bias_sum , pollen_biased_sum , pollen_tube_biased_sum )
        
        

        whole_genome_corr_subset$LG[whole_genome_corr_subset$LG == "LG4"] <- "XY"
        whole_genome_corr_subset$LG[whole_genome_corr_subset$LG == "LG1"] <- "A1"
        whole_genome_corr_subset$LG[whole_genome_corr_subset$LG == "LG3"] <- "A2"
        whole_genome_corr_subset$LG[whole_genome_corr_subset$LG == "LG5"] <- "A3"
        whole_genome_corr_subset$LG[whole_genome_corr_subset$LG == "LG2"] <- "A4"
        #  whole_genome_corr_subset$LG = factor(reshaped_map$LG, levels=c('A1','A2','A3','A4','XY'))
        
        
        
        
        chromosome_corr_subset<-subset(whole_genome_corr_subset, whole_genome_corr_subset$LG=="A1")
        main_title<-"A. A1 windowed partial correlations"
        
        chromosome_corr_subset<-subset(whole_genome_corr_subset, whole_genome_corr_subset$LG=="A2")
        main_title<-"B. A2 windowed partial correlations"
        
        chromosome_corr_subset<-subset(whole_genome_corr_subset, whole_genome_corr_subset$LG=="A3")
        main_title<-"C. A3 windowed partial correlations"
        
        chromosome_corr_subset<-subset(whole_genome_corr_subset, whole_genome_corr_subset$LG=="A4")
        main_title<-"D. A4 windowed partial correlations"
        
        chromosome_corr_subset<-subset(whole_genome_corr_subset, whole_genome_corr_subset$LG=="XY")
        main_title<-"E. XY windowed partial correlations"
        
        
        
        chromosome_corr_subset<-chromosome_corr_subset %>% dplyr::select(!LG)
        
        chromosome_partial_correlations_results<-data.frame(matrix(ncol = 7, nrow = 0))
        colnames(chromosome_partial_correlations_results)<-c("control","var_1","var_2","estimate","p","statistic","n")
        
        
        for (i in 2:length(colnames(chromosome_corr_subset))){
          #      print(i)
          # var_1<-i
          for (j in 2:length(colnames(chromosome_corr_subset))){
            if (i!=j){ 
              # var_2<-j
              #print(var_1)
              print(colnames(chromosome_corr_subset)[i])
              #print(var_2)
              print(colnames(chromosome_corr_subset)[j])
              # pcor.test(chromosome_corr_subset[,i],chromosome_corr_subset[,j],chromosome_corr_subset[,"n_genes_Yunchen"])
              df<-data.frame(chromosome_corr_subset[,i],chromosome_corr_subset[,j],chromosome_corr_subset[,"n_genes_Yunchen"])
              df<-df%>%drop_na(.)
              colnames(df)<-c(colnames(chromosome_corr_subset)[i],colnames(chromosome_corr_subset)[j],"n_genes_Yunchen")
              print( head(df))
              pcor_test_output<-pcor.test(df[,1], df[,2], df[3])
              print(pcor_test_output)
              outrow<-c(
                "n_genes_Yunchen",
                (colnames(chromosome_corr_subset)[i]),
                (colnames(chromosome_corr_subset)[j]),
                pcor_test_output$estimate,
                pcor_test_output$p.value,
                pcor_test_output$statistic,
                pcor_test_output$n
              )  
              outrow<-data.frame(t(unlist(outrow)))
              colnames(outrow)<-c("control","var_1","var_2","estimate","p","statistic","n")
              print(outrow)
              chromosome_partial_correlations_results<-rbind(chromosome_partial_correlations_results,outrow)
            }
          }
        }
        
        #View(chromosome_partial_correlations_results)
        
        
        
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^N_DNA_TEs_all", replacement = "DNA TEs per window")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^N_RNA_TEs_all", replacement = "RNA TEs per window")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^male_distortion_mean", replacement = "Male distortion")         
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^female_distortion_mean", replacement = "Female distortion")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^Mendelian_distortion_mean", replacement = "Mendelian distortion")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^crossovers_male_sum", replacement = "Male crossovers")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^crossovers_female_sum", replacement = "Female crossovers")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^crossovers_sex_averaged_sum", replacement = "Total crossovers")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^crossovers_diff", replacement = "Crossover difference")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^expressed_leaf_sum", replacement = "Leaf expressed genes")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^leaf_male_bias_sum", replacement = "Leaf male biased")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^leaf_female_bias_sum", replacement = "Leaf female biased")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^expressed_flower_sum", replacement = "Flower expressed genes")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^flower_male_bias_sum", replacement = "Flower male biased")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^flower_female_bias_sum", replacement = "Flower female biased")
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^pollen_biased_sum", replacement = "Pollen biased") 
        chromosome_partial_correlations_results[2:3] <- lapply(chromosome_partial_correlations_results[2:3], gsub, pattern = "^pollen_tube_biased_sum", replacement = "Pollen tube biased")
        
        chromosome_partial_correlations_results$estimate<-as.numeric(chromosome_partial_correlations_results$estimate)
        
        
        chrom_corr_plot<-chromosome_partial_correlations_results %>%
          mutate(var_1 = factor(var_1, levels=c(
            "DNA TEs per window", "RNA TEs per window", 
            "Male distortion", "Female distortion", "Mendelian distortion", "Leaf expressed genes", "Leaf male biased", "Leaf female biased", 
            "Flower expressed genes", "Flower male biased", "Flower female biased", 
            "Pollen biased", "Pollen tube biased"
          ))) %>%
          mutate(var_2 = factor(var_2, levels=c(
            "Male crossovers", "Female crossovers", "Total crossovers", "Crossover difference" #"Genes per window", "XO_diff",  "XO",  "XO_F",  "XO_M"
          ))) %>%  drop_na(., c(var_1,var_2))  %>%   ## to get the rect filled %>% 
          ggplot(aes(var_1,var_2, col=estimate)) +
          geom_tile(col="black", fill="white") + geom_point(aes(size = abs(estimate)), shape=15) + labs(x = NULL, y = NULL, col = "Pearson's\nCorrelation", title=paste(main_title)) + theme_classic() +
          scale_color_gradient2(mid="#FBFEF9",low="#0C6291",high="#A63446", limits=c(-1,1)) +
          scale_x_discrete(expand=c(0,0)) +
          scale_y_discrete(expand=c(0,0)) +
          #          scale_size(range=c(1,11), guide=NULL) + 
          scale_size(range=c(1,11), guide=NULL) + 
          theme(axis.text.x = element_text(angle=90))  
        
        chrom_corr_plot
        
        #ggsave("Figures/A1_partial_correlations.pdf", width = 8, height = 4, units = "in")
        #ggsave("Figures/A1_partial_correlations.png", width = 800, height = 400, units = "px", dpi=72)
        
        #ggsave("Figures/A2_partial_correlations.pdf", width = 8, height = 4, units = "in")
        #ggsave("Figures/A2_partial_correlations.png", width = 800, height = 400, units = "px", dpi=72)
        
        #ggsave("Figures/A3_partial_correlations.pdf", width = 8, height = 4, units = "in")
        #ggsave("Figures/A3_partial_correlations.png", width = 800, height = 400, units = "px", dpi=72)
        
        #ggsave("Figures/A4_partial_correlations.pdf", width = 8, height = 4, units = "in")
        #ggsave("Figures/A4_partial_correlations.png", width = 800, height = 400, units = "px", dpi=72)
        
        #ggsave("Figures/XY_partial_correlations.pdf", width = 8, height = 4, units = "in")
        #ggsave("Figures/XY_partial_correlations.png", width = 800, height = 400, units = "px", dpi=72)
        
        
        
        write.csv(chromosome_partial_correlations_results, "A1_windowed_partial_correlations.csv", quote=F, row.names=F)        
        write.csv(chromosome_partial_correlations_results, "A2_windowed_partial_correlations.csv", quote=F, row.names=F)        
        write.csv(chromosome_partial_correlations_results, "A3_windowed_partial_correlations.csv", quote=F, row.names=F)        
        write.csv(chromosome_partial_correlations_results, "A4_windowed_partial_correlations.csv", quote=F, row.names=F)        
        write.csv(chromosome_partial_correlations_results, "XY_windowed_partial_correlations.csv", quote=F, row.names=F)
        
        
        

    