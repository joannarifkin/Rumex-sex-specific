library(tidyverse)

snps<-read_table2("finalLG4SNPsAll", col_names = c("chr","pos","ref","alt","score","1","2","3","4","5","6","7","8","9","10","11","12"))
#theres a blank column at the end there, ignoring the error seems fine
xsnps <- read_table2("finalXrefxysnps", col_names = 
                       c("chr","pos","ref","alt","score","1","2","3","4","5","6","7","8","9","10","11","12")) %>%
  mutate(snp = "X")
ysnps <- read_table2("finalYrefxysnps", col_names = 
                       c("chr","pos","ref","alt","score","1","2","3","4","5","6","7","8","9","10","11","12")) 
  
#fixdifs <- bind_rows(xsnps,ysnps)

library(ggplot2)
####### Plot the densities of snps #####
ggplot(snps) + geom_histogram(aes(x=pos),binwidth=1e7) + 
  geom_histogram(data=xsnps, fill="green",aes(x=pos),binwidth=1e7) +
  geom_histogram(data=ysnps,fill="purple", aes(x=pos),binwidth=1e7) +
  ggtitle("snp snp snp") +
  xlab("Position in the genome") + 
  ylab("SNP density")

##### SNP proportions in 10Mb windows #####

# group by ranges
win_snp_num <- snps %>% mutate(win = cut(pos, breaks= seq(0, 240000000, by = 1e7, ordered_result=TRUE))) %>% 
  group_by(win) %>% 
  summarize(n=n())
win_xsnp <- xsnps %>% mutate(win = cut(pos, breaks= seq(0, 240000000, by = 1e7, ordered_result=T))) %>% 
  group_by(win) %>% 
  summarize(nx=n()) 
win_ysnp <- ysnps %>% mutate(win = cut(pos, breaks= seq(0, 240000000, by = 1e7, ordered_result=T))) %>% 
  group_by(win) %>% 
  summarize(ny=n()) 

fixdifs <-full_join(win_ysnp,win_xsnp, by = "win")
fixdifs[is.na(fixdifs)] = 0
allsnps <- full_join(win_snp_num, fixdifs,by = "win") %>% mutate(., propnx = nx/n ) %>% mutate(., propny = ny/n )
allsnps[is.na(allsnps)] = 0


ggplot(allsnps, aes(x=win))  + geom_line(color="green",aes(y=propnx,group=1)) + 
  geom_line(color="purple",aes(y=propny,group=1)) + 
  xlab("10 Mb window") + 
  ylab("SNP density")

##### SNP proportions in 1Mb windows #####

# group by ranges
win_snp_num_1mb <- snps %>% mutate(win = cut(pos, breaks= seq(0, 240000000, by = 1e6))) %>% 
  group_by(win) %>% 
  summarize(n=n())
win_xsnp_1mb <- xsnps %>% mutate(win = cut(pos, breaks= seq(0, 240000000, by = 1e6))) %>% 
  group_by(win) %>% 
  summarize(nx=n()) 
win_ysnp_1mb <- ysnps %>% mutate(win = cut(pos, breaks= seq(0, 240000000, by = 1e6, ordered_result	
=T))) %>% 
  group_by(win) %>% 
  summarize(ny=n()) 
fixdifs_1mb <-full_join(win_xsnp_1mb,win_ysnp_1mb, by = "win")
allsnps_1mb <- full_join(fixdifs_1mb, win_snp_num_1mb, by = "win") %>% 
  mutate(., propnx = nx/n ) %>% mutate(., propny = ny/n )



ggplot(allsnps_1mb, aes(x=win))  + geom_line(color="green",aes(y=propnx,group=1)) + 
  geom_line(color="purple",aes(y=propny,group=1)) + 
  xlab("10 Mb window") + 
  ylab("SNP density") 

##### fixed diffs irrespective of X or Y #####


allfixsnps <- allsnps %>% mutate(nfix = nx+ny) %>% mutate(., propnall = nfix/n ) %>% 
  rownames_to_column("win2") 
allfixsnps$win2 <- as.integer(allfixsnps$win2)
allfixsnps <- mutate(allfixsnps, Mb=win2*10)

# add dashed line to 212.1 MB
# add 
### line ###
ggplot(allfixsnps, aes(x=(Mb)))  + geom_line(aes(y=propnall, group=1)) +
  geom_vline(xintercept = 212.1, linetype= "dashed") +
  theme_bw(base_size = 14) +
  ggtitle("A") +
  labs(x="Position on sex chromosome (Mb)",y="Proportion fixed\nX-Y differences") +
  theme(strip.background =element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
ggsave("LG4_fixedSNPs_line.png", dpi=300, height=3.25, width=6.5, units="in")

#### lines colors ####
ggplot(allfixsnps, aes(x=(Mb)))  + 
  geom_line(group=1,color="darkgreen",aes(y=propnx)) +
  geom_line(group=1,color="purple",aes(y=propny)) +
  geom_vline(xintercept = 212.1, linetype= "dashed") + 
  theme_bw(base_size = 14) +
  ggtitle("B") +
  labs(x="Position on sex chromosome (Mb)",y="Proportion fixed\nX-Y differences") +
  theme(strip.background =element_rect(fill="white")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

  
ggsave("LG4_fixedSNPs_linesXY.png", dpi=300, height=3.25, width=6.5, units="in")

