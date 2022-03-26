install.packages("tidyverse")
install.packages("tableone")
install.packages("ggplot2")
library(plyr)
library(tidyverse) 
library(tableone)
library(multcomp)
library(MBESS)
library(flextable)
library(officer)
library(VennDiagram)
library(magrittr)
library(survival)
library(ggplot2)
source("~/Downloads/Ally data analysis/Scripts/makeDemo.R") 
source("~/Downloads/Ally data analysis/Scripts/runcontrast2.R")
source("~/Downloads/Ally data analysis/Scripts/returnAdj.R")
source("~/Downloads/Ally data analysis/Scripts/cohensd.R")
source('~/Downloads/Ally data analysis/Scripts/get_cor.R')

# Load the packages
cattable<-read.csv("~/Downloads/Ally data analysis/Data/CatTable.csv")
conttable<-read.csv("~/Downloads/Ally data analysis/Data/ContTable.csv")
resultstable<-read.csv("~/Downloads/Ally data analysis/All_Groups_CTQ.csv")
biotype<- read.csv("~/Downloads/Ally data analysis/Data/BSNIP2_Biotypes_rawdata.csv")


#substring SUBJID and remove duplicates for structural data
structure<-read.csv("~/Downloads/Ally data analysis/Data/allmeasures_05012020_brett.csv")
structure$SUBJID<- gsub("_.*", "", structure$SUBJID)
structure <- structure[!duplicated(structure$SUBJID),]
#read in CSV for all sites 
sites <- read.csv("~/Downloads/Ally data analysis/Data/allsites_05012020_brett.csv")


#clinical data cleaning and merging with structural data
bsnipclinical<- read.csv("~/Downloads/Ally data analysis/Data/bsnip2_ad_preliminary_20191001.csv", na.strings = c(999.0,999))
bsnip<-merge(x = bsnipclinical, y= structure, by.x="subject_id", by.y = "SUBJID", all = TRUE)
bsnip<- merge(x=bsnip, y= sites, by.x = "subject_id", by.y = "real_id", all = TRUE)

# removed 11  for missing LDPS score, 9 for missing age, 16 for Dx
bsnip <- bsnip[!is.na(bsnip$p1_andel_a) & !is.na(bsnip$age) & !bsnip$group=='' & !bsnip$group=='OTH',]
table(proband$q1_qual, proband$vh)
#calculated duration of illness
bsnip$doi<- bsnip$age - bsnip$szsadbponset

#created VH/AH, VH, AH, none, VH+, and VH- groups
bsnip$avh<- NA 
bsnip$avh[bsnip$p3_anhal_b=="1" &bsnip$p3_anhal_a=="1"] <- "AVH"
bsnip$avh[bsnip$p3_anhal_b=="1" &bsnip$p3_anhal_a=="0"] <- "VH"
bsnip$avh[bsnip$p3_anhal_b=="0" &bsnip$p3_anhal_a=="1"] <- "AH"
bsnip$avh[bsnip$p3_anhal_b=="0" &bsnip$p3_anhal_a=="0"] <- "nonAVH"
bsnip$avh[bsnip$group=="NC"] = "HC"

bsnip$vh<- NA
bsnip$vh[bsnip$p3_anhal_b == "1"] = "VH"
bsnip$vh[bsnip$p3_anhal_b == "0"&bsnip$group!="NC"] = "nonVH"
bsnip$vh[bsnip$group == "NC"] = "HC"
bsnip$vh2[bsnip$p3_anhal_b == "1"] = "VH"
bsnip$vh2[bsnip$p3_anhal_b == "0"&bsnip$group!="NC"] = "nonVH"

#separated race into AA, CA, and other 
bsnip$raz<-NA
bsnip$raz[bsnip$race=="AA"]<-"AA"
bsnip$raz[bsnip$race=="CA"]<-"CA"
bsnip$raz[bsnip$race== "AE" | bsnip$race == "AS"|bsnip$race == "MR"|bsnip$race == "OT"| bsnip$race == "NH"| bsnip$race== "UNK"] = "OTHER"

#separated subjects into HC and Probands with a total PANSS score
bsnip$phenotype<-NA
bsnip$phenotype[bsnip$group=="NC"]<- "HC"
bsnip$phenotype[bsnip$group!="NC"]<- "Probands"
bsnip$hall<- bsnip$p3_anhal_a + bsnip$p3_anhal_b + bsnip$p3_anhal_c + bsnip$p3_anhal_d + bsnip$p3_anhal_e
proband<- bsnip[bsnip$phenotype == "Probands",]

#visual hallucinator prevalence
prevalence <- 532/(532+632)*100  # 45.7%
prevalence2 <- 346/(346+492+286+40)*100  # 3.4% VH, 42.3% AVH, 24.6% AH, 29.7% nonAVH

bsnip2 <- bsnip[!is.na(bsnip$Left.Hippocampus),]
table(bsnip2$avh)

#create Table 1
listVars<-c("age", "sex", "raz", "group", "site", "sfs_total",  "GAF","ConsensusWRAT","sub_edscale", "pt_hollscore","panss_gentotal", "panss_postotal", "panss_negtotal", "madrs_total", "young_total", "p3_anhal_a", "p3_anhal_c","p3_anhal_d", "p3_anhal_e", "p1_andel_a", "p1_andel_b", "p1_andel_c", "p1_andel_d", "p1_andel_e", "p1_andel_f", "med_antip_any_pri", "med_antip_firstgen_sec", "med_antip_secondgen_sec","doi", "cas_total","hall", "med_avgdailycpz")
catVars <- c("sex", "raz", "group", "site", "educ","hollclass","p3_anhal_a", "p3_anhal_c","p3_anhal_d", "p3_anhal_e", "p1_andel_a", "p1_andel_b", "p1_andel_c", "p1_andel_d", "p1_andel_e", "p1_andel_f","med_antip_any_pri", "med_antip_firstgen_sec", "med_antip_secondgen_sec", "Biotype")
ctqVars <- c('ctq_total','ctqscore_val','ctqscore_ea','ctqscore_en','ctqscore_pa','ctqscore_pn','ctqscore_sa')
listVars2 <- c("age", "sex", "raz", "quantile", "Left_V1_exvivo_area", "Left_V1_exvivo_thickness", "Left_V1_exvivo_volume", 
               "Left_V2_exvivo_area", "Left_V2_exvivo_thickness", "Left_V2_exvivo_volume",
               "Right_V1_exvivo_area", "Right_V1_exvivo_thickness", "Right_V1_exvivo_volume",
               "Right_V2_exvivo_area", "Right_V2_exvivo_thickness", "Right_V2_exvivo_volume",
               "Left_MT_exvivo_area", "Left_MT_exvivo_thickness", "Left_MT_exvivo_volume",
               "Right_MT_exvivo_area", "Right_MT_exvivo_thickness", "Right_MT_exvivo_volume",
               "Left_fusiform_area", "Left_fusiform_thickness", "Left_fusiform_volume",
               "Right_fusiform_area", "Right_fusiform_thickness", "Right_fusiform_volume",
               "Left_lingual_area", "Left_lingual_thickness", "Left_lingual_volume",
               "Right_lingual_area", "Right_lingual_thickness", "Right_lingual_volume",
               "Left.Thalamus", "Right.Thalamus", "Left.Hippocampus", "Right.Hippocampus")
catVars2 <- c("sex", "raz")
table1 <-CreateTableOne(vars=listVars, data =proband, factorVars = catVars, strata = c("vh"))
hctable <-CreateTableOne(vars=listVars, data =bsnip, factorVars = catVars, strata = c("vh"))
avhtable<-CreateTableOne(vars=listVars, data =proband, factorVars = catVars, strata = c("avh"))
ctqtable<-CreateTableOne(vars=ctqVars, data =proband, factorVars = catVars, strata = c("avh"))
strtrable<-CreateTableOne(vars=listVars2, data =bsnip, factorVars = catVars2, strata = c("avh"))
prostrtrable<-CreateTableOne(vars=listVars2, data =proband, factorVars = catVars2, strata = c("avh"))

table(bsnip2$quantile, bsnip2$phenotype)

#convert table 1 to Excel
print(table1,  quote = TRUE, noSpaces = TRUE)
print(hctable,  quote = TRUE, noSpaces = TRUE)
print(avhtable,  quote = TRUE, noSpaces = TRUE)
print(ctqtable,  quote = TRUE, noSpaces = TRUE)

# Stepwise Regression
library(MASS)
fit <- lm(ctq_total ~   sex + race + site + med_antip_secondgen_sec + pt_hollscore, data=proband)
step <- stepAIC(fit, direction="both")

# bsnip_proband_bacs <- proband[(is.na(proband$pt_hollclass)==FALSE),]
# fit <- lm(proband$BACS_COMP_z ~   pt_hollclass+ race+ site + med_antip_secondgen_sec,data=proband)
# step <- stepAIC(fit, direction="both")

cattable$Variable<- as.character(cattable$Variable)

label.cat3 <- data.frame(Variable = c("Auditory", "Tactile", "Olfactory", "Gustatory", "Paranoid", "Grandiose", "Somatic", "Religious"),
                       Percentage = c(95, 48, 30, 14, 90, 51, 31, 40), Group = "VH")
label.cat2 <-data.frame(Variable = c("Other"), Percentage = c(30), Group = "VH")

figure1<- ggplot(data = cattable, aes(x=Variable, y=Percentage, fill = Group)) + geom_bar(stat="identity", position=position_dodge()) + 
  # scale_fill_manual(values=c("red", "blue"))+
  ggtitle('Comparing Clinical Variables Between non-VH and VH Groups') +
  geom_text(data = label.cat3, label = "***", position=position_dodge(.9), size = 6) + 
  geom_text(data = label.cat2, label = "**", position=position_dodge(.9), size = 6)+
  xlab("")+theme(legend.text= element_text(size=15), axis.text.x=element_text(size = 15),
  axis.text.y=element_text(size = 16, hjust = 1.0), axis.title.y =element_text(size=24),
  title=element_text(size=24), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))

figure1gray<- ggplot(data = cattable, aes(x=Variable, y=Percentage, fill = Group)) + geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_grey()+
  ggtitle('Comparing Clinical Variables Between non-VH and VH Groups') +
  geom_text(data = label.cat3, label = "***", position=position_dodge(.9), size = 6) + 
  geom_text(data = label.cat2, label = "**", position=position_dodge(.9), size = 6)+
  xlab("")+theme(legend.text= element_text(size=15), axis.text.x=element_text(size = 15),
                 axis.text.y=element_text(size = 16, hjust = 1.0), axis.title.y =element_text(size=24),
                 title=element_text(size=24), plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
table(proband$group, proband$vh2)

#Rename factor levels of group
levels(proband$group)<-c("","BPP","NC","OTH","SAD","SZ")
proband$group<-factor(proband$group) #drop factor levels

#add y- axis label and numbers in boxes
#Change BT names
proband$biotype2<-proband$Biotype #create new cleaned version of biotype
proband$biotype2[proband$Biotype=="1"]<-"BT1" #re-code erroneous biotype values
proband$biotype2[proband$Biotype=="2"]<-"BT2"
proband$biotype2[proband$Biotype=="5"]<-"BT3"
table(proband$biotype2) #verify re-coding
table(proband$biotype2, proband$Biotype)#check re-coding of biotype
proband$biotype2<-factor(proband$biotype2) #drop as a factor level
#Biotype figure
biotype.vh.tab<-table(proband$vh2,proband$biotype2)
biotype.vh.tab
# Create barplot with count
barplot(biotype.vh.tab, beside = F, legend.text = T)
#Diagnosis figure
group.vh.tab<-table(proband$p3_anhal_b, proband$group)
group.vh.tab
# Create barplot with count
barplot(group.vh.tab, beside = F, legend.text = T)
combinedplot<- cbind(biotype.vh.tab,group.vh.tab)
row<-c(92, 85)
row1<- c(92, 85, 114, 178, 184, 270)
row2 <- c(112, 111, 125, 81, 242, 209)
plot <- barplot(combinedplot,beside=F, legend.text=T, main = "Comparing Biotypes and Diagnoses between VH and non-VH",
        ylab = "Count", 
        xlim=c(0, ncol(combinedplot) + 2), ylim = c(0, 500),
        args.legend=list(
          x=ncol(combinedplot) + 3,
          y=max(colSums(combinedplot)),
          bty = "n"))
text(x = plot, y = combinedplot, pos = 3, cex = 0.8, col = "black")
text(plot, 5, round(row, 1),cex=1,pos=3) 
text(plot, 5, round(row1, 1),cex=1,pos=3) 
text(plot, 150, round(row2, 1), cex=1, pos=3)

library(plyr)
library("ggpubr")
# Create the barplot
proband$Biotype <- as.factor(proband$Biotype)
proband$Biotype<-factor(proband$Biotype)
btplot <-ggplot(data=subset(proband, !is.na(Biotype)), aes(x=Biotype, fill =  vh2 )) +
  labs(fill = "Group")+
  geom_bar(stat="count", width = 0.6) +theme(text = element_text(size=16))+ geom_text(size =6, aes(label=..count..),stat="count",position=position_stack(0.5)) +
  scale_x_discrete(breaks=c(1, 2,3,4,5),labels=c("BT1", "BT2", "BT3", "BT4", "BT3"))
table(proband$Biotype)
groupplot<-ggplot(data=proband, aes(x=group, fill = vh2 )) +
  labs(fill = "Group")+theme(text = element_text(size=16))+
  geom_bar(stat="count", width = 0.6)  + scale_x_discrete(name = "Diagnosis",labels=c("BPP", "SAD", "SZ"))+geom_text(size = 6, aes(label=..count..),stat="count",
                                                                                  position=position_stack(0.5))
btgroupplot<- rbind(btplot,groupplot)
combined<- ggplot(combinedplot, aes(x=ISSUE_DATE, group=group, col=group, fill=group)) +
  geom_bar(stat="count",position = "dodge")

#display for panss
proband$vh <- factor(proband$vh, levels = c('nonVH', 'VH'))
pheno.mat2<- rbind(c(1,-1))
rownames(pheno.mat2)<- c("VH-nonVH")
covars_r<-c('age','sex','raz','site')
symp <- c('panss_postotal','panss_gentotal','panss_negtotal','cas_total','young_total','madrs_total')
results_panss<- RunMat(symp, pheno.mat2, proband, 'vh', covars_r)
results_panss$p_adj<-p.adjust(p=results_panss$p, method = 'fdr')
results_panss<-get_runcontrast_d(proband,results_panss, 'vh', covars_r)
results_panss[results_panss$p_adj<0.01,]
#adding CPZ as covars
covars_cpz<-c('age','sex','raz','site', 'med_avgdailycpz')
panss<-names(proband)[grep("panss_pos|panss_gen|panss_tot", names(proband))]
symp <- c('panss_postotal','panss_gentotal','panss_negtotal','cas_total','young_total','madrs_total')
results_pansscpz<- RunMat(symp, pheno.mat2, proband, 'vh', covars_cpz)
results_pansscpz$p_adj<-p.adjust(p=results_pansscpz$p, method = 'fdr')
results_pansscpz<-get_runcontrast_d(proband,results_pansscpz, 'vh', covars_cpz)
results_pansscpz[results_pansscpz$p_adj<0.9,]
#adding q1_qual as covars
covars_qual<-c('age','sex','raz','site', 'q1_qual')
panss<-names(proband)[grep("panss_pos|panss_gen|panss_tot", names(proband))]
symp <- c('panss_postotal','panss_gentotal','panss_negtotal','cas_total','young_total','madrs_total')
results_panssqual<- RunMat(symp, pheno.mat2, proband, 'vh', covars_qual)
results_panssqual$p_adj<-p.adjust(p=results_panssqual$p, method = 'fdr')
results_panssqual<-get_runcontrast_d(proband,results_panssqual, 'vh', covars_qual)
results_panssqual[results_panssqual$p_adj<0.9,]

fp <- ggplot(data=results_panss, aes(x=region, y=d, ymin=CI_lower, ymax=CI_upper, colour=contr)) 
fp+geom_pointrange(size=1.5, position = position_dodge(width=0.5)) + scale_colour_manual(values=c('blue'))+
geom_hline(yintercept=0, lty=2, size=1.5, color='red')+  # add a dotted line at x=1 after flip
coord_flip() + ylim(c(-0.15,0.5))+ # flip coordinates (puts labels on y axis)
xlab("") + ylab("Cohen's d (95% CI)") + labs(color = "Contrast")+
theme_bw()+ggtitle('Comparing Clinical Variables Between VH and non-VH')+
theme(legend.text=element_text(size=16), legend.title = element_text(size=18),
axis.text.x=element_text(size = 18),axis.text.y=element_text(size = 18),
axis.title.y = element_text(size=22),axis.title.x = element_text(size=22),
plot.title = element_text(size=27, hjust=0.5)) + 
scale_x_discrete(labels=c("madrs_total"="MADRS","young_total"='YMRS',"cas_total"="CAS", "panss_negtotal"= "PANSS Negative", "panss_postotal"="PANSS Positive","panss_gentotal"='PANSS General'))     
        
#run contrast for cognition
cog<-names(bsnip)[grep("BACS", names(bsnip))]
covars_cog<-c('raz','site')
results_cog<- RunMat(cog, pheno.mat2, proband, 'vh', covars_cog)
results_cog$p_adj<-p.adjust(p=results_cog$p, method = 'fdr')
results_cog<-get_runcontrast_d(proband,results_cog, 'vh', covars_cog)
results_cog[results_cog$p_adj<0.1,]

# group constrasts for refined analysis   
bsnip$avh <- factor(bsnip$avh, levels = c('VH', 'AH', 'AVH','nonAVH'))
pheno.mat<- rbind(c(1,-1,0,0), c(1,0,-1,0), c(1,0,0,-1),c(0,1,-1,0),c(0,1,0,-1),c(0,0,1,-1))
rownames(pheno.mat)<- c("VH-AH", "VH-AVH", "VH-nonAVH","AH-AVH","AH-nonAVH","AVH-nonAVH")
covars_r<-c('age','sex','raz','site')
results_cog<- RunMat(symp, pheno.mat, bsnip, 'avh', covars_r)
results_cog$p_adj<-p.adjust(p=results_cog$p, method = 'fdr')
results_cog<-get_runcontrast_d(bsnip,results_cog, 'avh', covars_r)
results_cog[results_cog$p_adj<0.05,]

bsnip$avh <- factor(bsnip$avh, levels = c('VH', 'AH', 'AVH','nonAVH'))
pheno.mat<- rbind(c(1,-1,0,0), c(1,0,-1,0), c(1,0,0,-1),c(0,1,-1,0),c(0,1,0,-1),c(0,0,1,-1))
rownames(pheno.mat)<- c("VH-AH", "VH-AVH", "VH-nonAVH","AH-AVH","AH-nonAVH","AVH-nonAVH")
covars_r<-c('age','sex','raz','site')
results_cog<- RunMat(symp, pheno.mat, bsnip, 'avh', covars_r)
results_cog$p_adj<-p.adjust(p=results_cog$p, method = 'fdr')
results_cog<-get_runcontrast_d(bsnip,results_cog, 'avh', covars_r)
results_cog[results_cog$p_adj<0.05,]

# fine grain comparison for cognition
results_cog2<- RunMat(cog, pheno.mat, bsnip, 'avh', covars_cog)
results_cog2$p_adj<-p.adjust(p=results_cog2$p, method = 'fdr')
results_cog2<-get_runcontrast_d(bsnip,results_cog2, 'avh', covars_cog)
results_cog2[results_cog2$p_adj<0.1,]
 
ctq<- c("ctqscore_ea","ctqscore_pa","ctqscore_sa","ctqscore_en","ctqscore_pn","ctqscore_val","ctq_total")
bacs3_corr<- get_correlations(bsnip, meas1 = "BACS_COMP_z", meas2 = ctq, covars1 = covars_cog, covars2 = covars_r, method = 'spearman')
bacs3_corr[bacs3_corr$adjusted_p_value<0.05,]

bacs3_corr_pro<- get_correlations(proband, meas1 = "BACS_COMP_z", meas2 = ctq, covars1 = covars_cog, covars2 = covars_r, method = 'spearman')
bacs3_corr_pro[bacs3_corr_pro$adjusted_p_value<0.05,]

bacs3_corr_hc<- get_correlations(hc, meas1 = "BACS_COMP_z", meas2 = ctq, covars1 = covars_cog, covars2 = covars_r, method = 'spearman')
bacs3_corr_hc[bacs3_corr_hc$adjusted_p_value<0.05,]

##### Paulo stopped here   ##### 
##### ##### ##### ##### ##### ##### 
##### ##### ##### ##### ##### ##### 

#run contrast
# bsnip$vh <- factor(bsnip$vh, levels = c('HC', 'VH', 'nonVH'))
# pheno.mat<- rbind(c(-1,1,0), c(-1,0,1), c(0,-1,1))
# rownames(pheno.mat)<- c("HC-nonVH", "HC-VH", "nonVH-VH")
# covars_r<-c('raz','site')
# cog<-names(bsnip)[grep("BACS", names(bsnip))]
# # cog<-c("BACS_COMP_z","GAF",'sfs_total','ctq_total')
# results_cog<- RunMat(cog, pheno.mat, bsnip, 'vh', covars_r)
# results_cog$p_adj<-p.adjust(p=results_cog$p, method = 'fdr')
# results_cog<-get_runcontrast_d(bsnip,results_cog, 'vh', covars_r)
# results_cog[results_cog$p_adj<0.05,]
#       
# bsnip$avh <- factor(bsnip$avh, levels = c('HC', 'VH', 'AH', 'AVH','nonAVH'))
# pheno.mat<- rbind(c(-1,1,0,0,0), c(-1,0,1,0,0), c(-1,0,0,1,0), c(-1,0,0,0,1))
# rownames(pheno.mat)<- c("HC-VH", "HC-AH", "HC-AVH", "HC-nonAVH")
# covars_r<-c('raz','site')
# cog<-names(bsnip)[grep("BACS", names(bsnip))]
# # cog<-c("BACS_COMP_z","GAF",'sfs_total','ctq_total')
# results_cog<- RunMat(cog, pheno.mat, bsnip, 'avh', covars_r)
# results_cog$p_adj<-p.adjust(p=results_cog$p, method = 'fdr')
# results_cog<-get_runcontrast_d(bsnip,results_cog, 'avh', covars_r)
# results_cog[results_cog$p_adj<0.05,]
        
#boxplot
# ggplot(bsnip, aes(x = avh, y = BACS_COMP_z, group=avh))+geom_boxplot(aes(fill=vh))+geom_jitter(width=0.3, size=3)+
# scale_fill_brewer(palette="Blues")+stat_summary(fun.y = "mean", colour = "red", size = 5, geom = "point")+
# ggtitle('Analytes in Diagnostic Groups & Controls')+
# theme(legend.position="none",axis.text.x=element_text(size = 26),
# axis.text.y=element_text(size = 22), axis.title.y =element_text(size=28),
# title=element_text(size=24), plot.title = element_text(hjust = 0.5)) + xlab('')
        

proband$avh <- factor(proband$avh, levels = c('VH', 'AH', 'VH/AH'))
pheno.mat2<- rbind(c(-1,1,0), c(-1,0,1), c(0,-1,1))
rownames(pheno.mat2)<- c("VH-AH", "VH-VH/AH", "AH-VH/AH")
covars_r<-c('age','sex','raz','site')
panss<-names(proband)[grep("panss_pos|panss_gen|panss_tot", names(proband))]
results_panss<- RunMat(panss, pheno.mat2, proband, 'avh', covars_r)
results_panss$p_adj<-p.adjust(p=results_panss$p, method = 'fdr')
results_panss<-get_runcontrast_d(proband,results_panss, 'avh', covars_r)
results_panss[results_panss$p_adj<0.05,]

#display for ctq
bsnip$vh <- factor(bsnip$vh, levels = c('HC', 'VH', 'nonVH'))
pheno.mat<- rbind(c(-1,1,0), c(-1,0,1), c(0,-1,1))
rownames(pheno.mat)<- c("HC-nonVH", "HC-VH", "nonVH-VH")
covars_r<-c('raz','site')
ctq<-names(bsnip)[grep("ctqscore_e|ctqscore_s|ctqscore_p|ctq_total", names(bsnip))]
results_ctq<- RunMat(ctq, pheno.mat, bsnip, 'vh', covars_r)
results_ctq$p_adj<-p.adjust(p=results_ctq$p, method = 'fdr')
results_ctq<-get_runcontrast_d(bsnip,results_ctq, 'vh', covars_r)
results_ctq[results_ctq$p_adj<0.01,]

ctq<-names(bsnip)[grep("ctqscore_e|ctqscore_s|ctqscore_p|ctq_total", names(bsnip))]
results_ctq<- RunMat(ctq, pheno.mat, bsnip, 'avh', covars_r)
results_ctq$p_adj<-p.adjust(p=results_ctq$p, method = 'fdr')
results_ctq<-get_runcontrast_d(bsnip,results_ctq, 'avh', covars_r)
results_ctq[results_ctq$p_adj<0.05,]

gaf<-names(bsnip)[grep("GAF", names(bsnip))]
results_gaf<- RunMat(ctq, pheno.mat, bsnip, 'vh', covars_r)
results_gaf$p_adj<-p.adjust(p=results_gaf$p, method = 'fdr')
results_gaf<-get_runcontrast_d(bsnip,results_gaf, 'vh', covars_r)
results_gaf[results_gaf$p_adj<0.05,]

HC<- bsnip[bsnip$group=="NC",]
VH<- proband[proband$p3_anhal_b == "1",]
nonVH<- proband[!proband$hall=="0"&proband$p3_anhal_b=="0",]

cog<- c('BACS_COMP_z', 'BACS_Verb_Mem_z', 'BACS_Verb_Flu_z')
hall<- c('p3_duration', 'p3_severity')
covars_cog<-c('raz','site')
covars_hall<-c('age', 'sex','raz', 'site')
hall_corr<-get_correlations(df = VH, meas1 = cog, meas2 = hall, covars1 = covars_cog, covars2 = covars_hall, method = 'spearman')
hall_corr2<-get_correlations(df = nonVH, meas1 = cog, meas2 = hall, covars1 = covars_cog, covars2 = covars_hall, method = 'spearman')

ctq<- c('ctq_total', 'ctqscore_ea', 'ctqscore_sa', 'ctqscore_pa', 'ctqscore_en', 'ctqscore_pn')
hall<- c('hall')
covars_ctq<-c('age', 'sex', 'raz','site')
hall_ctq<-get_correlations(df = VH, meas1 = ctq, meas2 = hall, covars1 = covars_ctq, covars2 = covars_hall, method = 'spearman')
hall_ctq2<-get_correlations(df = nonVH, meas1 = ctq, meas2 = hall, covars1 = covars_ctq, covars2 = covars_hall, method = 'spearman')

ctq<- c('ctq_total', 'ctqscore_ea', 'ctqscore_sa', 'ctqscore_pa', 'ctqscore_en', 'ctqscore_pn')
hall<- c('hall')
covars_ctq<-c('age', 'sex', 'raz','site')
hall_ctq<-get_correlations(df = AH, meas1 = ctq, meas2 = hall, covars1 = covars_ctq, covars2 = covars_hall, method = 'spearman')
hall_ctq2<-get_correlations(df = nonAH, meas1 = ctq, meas2 = hall, covars1 = covars_ctq, covars2 = covars_hall, method = 'spearman')

ctq<- c('ctq_total', 'ctqscore_ea', 'ctqscore_sa', 'ctqscore_pa', 'ctqscore_en', 'ctqscore_pn')
hall<- c('hall')
covars_ctq<-c('age', 'sex', 'raz','site')
hall_ctq<-get_correlations(df =proband, meas1 = ctq, meas2 = hall, covars1 = covars_ctq, covars2 = covars_hall, method = 'spearman')
#hall_ctq2<-get_correlations(df = nonVH, meas1 = ctq, meas2 = hall, covars1 = covars_ctq, covars2 = covars_hall, method = 'spearman')

cog<- c('BACS_COMP_z', 'BACS_Verb_Mem_z', 'BACS_Verb_Flu_z', 'BACS_Dig_Seq_z', 'BACS_Sym_Cod_z', 'BACS_Tok_Mot_z', 'BACS_Tower_z')
ctq<- c('ctq_total', 'ctqscore_ea', 'ctqscore_sa', 'ctqscore_pa', 'ctqscore_en', 'ctqscore_pn')
covars_cog<-c('raz','site')
bacs_corr<-get_correlations(df = proband, meas1 = cog, meas2 = ctq, covars1 = covars_cog, covars2 = covars_ctq, method = 'spearman')bacs_corr2<-get_correlations(df = nonVH, meas1 = cog, meas2 = ctq, covars1 = covars_cog, covars2 = covars_ctq, method = 'spearman')
bacs_corr3<-get_correlations(df = HC, meas1 = cog, meas2 = ctq, covars1 = covars_cog, covars2 = covars_ctq, method = 'spearman')
bacs_corr4<-get_correlations(df = VH, meas1 = cog, meas2 = ctq, covars1 = covars_cog, covars2 = covars_ctq, method = 'spearman')

cog<- c('BACS_COMP_z')
ctq<- c('ctq_total', 'ctqscore_ea', 'ctqscore_sa', 'ctqscore_pa', 'ctqscore_en', 'ctqscore_pn')
covars_cog<-c('raz','site')
bacs2_corr<- get_correlations(df = proband, meas1 = cog, meas2 = ctq, covars1 = covars_cog, covars2 = covars_ctq, method = 'spearman')
bacs3_corr<- get_correlations(df = HC, meas1 = cog, meas2 = ctq, covars1 = covars_cog, covars2 = covars_ctq, method = 'spearman')
bacs4_corr<- get_correlations(df = bsnip, meas1 = cog, meas2 = ctq, covars1 = covars_cog, covars2 = covars_ctq, method = 'spearman')
bacs4_corr[bacs4_corr$adjusted_p_value<0.001,]

cog<- c('BACS_COMP_z')
hall<- c('hall')
covars_cog<-c('raz','site')
bacs3_corr<- get_correlations(df = proband, meas1 = cog, meas2 = hall, covars1 = covars_cog, covars2 = covars_ctq, method = 'spearman')

proband$BACS_Verb_Mem_z_adj<- returnAdj(proband, 'BACS_Verb_Mem_z', covars = covars_cog)
proband$p3_duration_adj<- returnAdj(proband, 'p3_duration', covars = covars_r)
proband$hall_adj<- returnAdj(proband, 'hall', covars = covars_ctq)
p1<-ggplot(proband,aes(p3_duration_adj, BACS_Verb_Mem_z_adj, color = vh)) +geom_point(aes(colour=vh),size=4) + 
  geom_smooth(method="lm",size=1.5,se=F) + labs(color = "Group") +
  xlab("Adjusted Hallucination Duration") + ylab("Adjusted BACS Verbal Memory") +
  theme_linedraw() +
  theme(legend.text=element_text(size=25), legend.title = element_text(size=30),
        axis.text.x=element_text(size = 20),axis.text.y=element_text(size = 20),
        axis.title.y = element_text(size=25),axis.title.x = element_text(size=25))
p1

proband$ctq_total_adj<- returnAdj(proband, 'ctq_total', covars = covars_ctq)
p2<-ggplot(proband,aes(p3_duration_adj, ctq_total_adj, color = vh)) +geom_point(aes(colour=vh),size=4) + 
  geom_smooth(method="lm",size=1.5,se=F) +
  xlab("Adjusted Hallucination Duration") + ylab("Adjusted CTQ Total") +
  theme_linedraw() + labs(color = "Group") +
  theme(legend.text=element_text(size=25), legend.title = element_text(size=30),
        axis.text.x=element_text(size = 20),axis.text.y=element_text(size = 20),
        axis.title.y = element_text(size=25),axis.title.x = element_text(size=25))
p2  

proband$ctqscore_sa_adj<- returnAdj(proband, 'ctqscore_sa', covars = covars_ctq)
p3<-ggplot(proband,aes(hall_adj, ctqscore_sa_adj, color = vh)) +geom_point(aes(colour=vh),size=4) + 
  geom_smooth(method="lm",size=1.5,se=F) +
  xlab("Adjusted Hallucination Total") + ylab("Adjusted CTQ Total") +
  theme_linedraw() + labs(color = "Group") +
  theme(legend.text=element_text(size=25), legend.title = element_text(size=30),
        axis.text.x=element_text(size = 20),axis.text.y=element_text(size = 20),
        axis.title.y = element_text(size=25),axis.title.x = element_text(size=25))

figure1<- ggplot(data = cattable, aes(x=Variable, y=Percentage, fill = Group)) + geom_bar(stat="identity", position=position_dodge()) + 
  ggtitle('Comparing Clinical Variables Between non-VH and VH Groups') +
  geom_text(data = label.cat3, label = "***", position=position_dodge(.9)) + 
  geom_text(data = label.cat2, label = "**", position=position_dodge(.9))+
  xlab("")+theme(legend.text= element_text(size=15), axis.text.x=element_text(size = 15),
                 axis.text.y=element_text(size = 16, hjust = 1.0), axis.title.y =element_text(size=24),
                 title=element_text(size=24), plot.title = element_text(hjust = 0.5))

data <- bsnip[bsnip$p3_anhal_a == "1" & bsnip$p3_anhal_b == "1" & bsnip$p3_anhal_d == "1" & bsnip$p3_anhal_e=="1",]

##########################structural


#density plot
hist(proband$ctq_total)
hist(bsnip$ctq_total)
hist(hc$ctq_total)
summary(proband$ctq_total)
proband$ctq_total
bsnip$q1
hc$ctq_total
hc$quantile <- NA
bsnip <- bsnip[is.na(bsnip$ctq_total)==FALSE,]
#used summary function for quantiles. used quantiles generated from first study and applied to HC
bsnip$quantile[bsnip$ctq_total<=45] = "q1"
bsnip$quantile[bsnip$ctq_total >=46 & bsnip$ctq_total<=56] = "q2"
bsnip$quantile[bsnip$ctq_total >=57 & bsnip$ctq_total<=71] = "q3"
bsnip$quantile[bsnip$ctq_total >72] = "q4"

#bsnip
table(bsnip$quantile, bsnip$phenotype)
summary(bsnip$quantile)
tapply(bsnip$ctq_total, bsnip$quantile, mean, na.rm=T)

#comparing min and max for ctq total 
proband <- bsnip[bsnip$phenotype == "Probands",]
hc<- bsnip[bsnip$phenotype == "HC",]
tapply(proband$ctq_total, proband$quantile, mean, na.rm=T)
tapply(hc$ctq_total, hc$quantile, mean, na.rm=T)
summary(hc$ctq_total)
summary(proband$ctq_total)


bsnip_pro <- bsnip_main[!(bsnip_main$HAL_TYPE=='2'),]

listVars2 <- c("age", "sex", "raz", "quantile", "Left_V1_exvivo_area", "Left_V1_exvivo_thickness", "Left_V1_exvivo_volume", 
               "Left_V2_exvivo_area", "Left_V2_exvivo_thickness", "Left_V2_exvivo_volume",
               "Right_V1_exvivo_area", "Right_V1_exvivo_thickness", "Right_V1_exvivo_volume",
               "Right_V2_exvivo_area", "Right_V2_exvivo_thickness", "Right_V2_exvivo_volume",
               "Left_MT_exvivo_area", "Left_MT_exvivo_thickness", "Left_MT_exvivo_volume",
               "Right_MT_exvivo_area", "Right_MT_exvivo_thickness", "Right_MT_exvivo_volume",
               "Left_fusiform_area", "Left_fusiform_thickness", "Left_fusiform_volume",
               "Right_fusiform_area", "Right_fusiform_thickness", "Right_fusiform_volume",
               "Left_lingual_area", "Left_lingual_thickness", "Left_lingual_volume",
               "Right_lingual_area", "Right_lingual_thickness", "Right_lingual_volume",
               "Left.Thalamus", "Right.Thalamus", "Left.Hippocampus", "Right.Hippocampus")
catVars2 <- c("sex", "raz")
avhstrtrable<-CreateTableOne(vars=listVars2, data =bsnip, factorVars = catVars2, strata = c("avh"))
vhstrtrable<-CreateTableOne(vars=listVars2, data =bsnip, factorVars = catVars2, strata = c("vh"))
avhprostrtrable<-CreateTableOne(vars=listVars2, data =proband, factorVars = catVars2, strata = c("avh"))
vhprostrtrable<-CreateTableOne(vars=listVars2, data =proband, factorVars = catVars2, strata = c("vh"))

proband <- proband[is.na(proband$ctq_total)==FALSE,]
proband$CTQ_QUANTILES <-  cut(proband$ctq_total, breaks=c(quantile(proband$ctq_total, probs = seq(0, 1, by = 0.25))), 
                                         labels=c("Q1","Q2","Q3","Q4"), include.lowest = TRUE, na.rm=F)
table(bsnip$quantile, exclude = NULL)
bsnip$quantile <- as.numeric(as.factor(bsnip$quantile))
table(bsnip$quantile)
class(bsnip$quantile)

q1 <- proband[proband$CTQ_QUANTILES == "Q1",]
summary(q1$ctq_total)
q2 <- proband[proband$CTQ_QUANTILES == "Q2",]
summary(q2$ctq_total)

table(proband$CTQ_QUANTILES)
summary(proband$CTQ_QUANTILES)  

table(bsnip$vh)
#run contrast for brain structures with VH and area
bsnip$vh <- factor(bsnip$vh, levels = c('HC', 'VH', 'nonVH'))
pheno.mat<- rbind(c(-1,1,0), c(-1,0,1), c(0,-1,1))
rownames(pheno.mat)<- c("HC-nonVH", "HC-VH", "nonVH-VH")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
str<-names(bsnip)[grep("Left_V1_exvivo_area|Left_V2_exvivo_area|Right_V1_exvivo_area|Right_V2_exvivo_area
|Left_MT_exvivo_area|Right_MT_exvivo_area|Left_fusiform|Right_fusiform_area|Left_lingual_area|Right_lingual_area|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'vh', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'vh', covars_r)
results_str[results_str$p_adj<0.1,]

#run contrast for brain structures with VH and volume
bsnip$vh <- factor(bsnip$vh, levels = c('HC', 'VH', 'nonVH'))
pheno.mat<- rbind(c(-1,1,0), c(-1,0,1), c(0,-1,1))
rownames(pheno.mat)<- c("HC-nonVH", "HC-VH", "nonVH-VH")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
str<-names(bsnip)[grep("Left_V1_exvivo_volume|Left_V2_exvivo_volume|Right_V1_exvivo_volume|Right_V2_exvivo_volume
|Left_MT_exvivo_volume|Right_MT_exvivo_volume|Left_fusiform_exvivo_volume|Right_fusiform_volume|Left_lingual_volume|Right_lingual_volume|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'vh', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'vh', covars_r)
results_str[results_str$p_adj<0.1,]

#run contrast for brain structures with VH and thickness
bsnip$vh <- factor(bsnip$vh, levels = c('HC', 'VH', 'nonVH'))
pheno.mat<- rbind(c(-1,1,0), c(-1,0,1), c(0,-1,1))
rownames(pheno.mat)<- c("HC-nonVH", "HC-VH", "nonVH-VH")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
str<-names(bsnip)[grep("Left_V1_exvivo_thickness|Left_V2_exvivo_thickness|Right_V1_exvivo_thickness|Right_V2_exvivo_thickness
|Left_MT_exvivo_thickness|Right_MT_exvivo_thickness|Left_fusiform_thickness|Right_fusiform_thickness|Left_lingual_thickness|Right_lingual_thickness|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'vh', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'vh', covars_r)
results_str[results_str$p_adj<0.1,]

#run contrast with avh
bsnip$avh <- factor(bsnip$avh, levels = c('VH', 'AH', 'AVH','nonAVH'))
pheno.mat<- rbind(c(1,-1,0,0), c(1,0,-1,0), c(1,0,0,-1),c(0,1,-1,0),c(0,1,0,-1),c(0,0,1,-1))
rownames(pheno.mat)<- c("VH-AH", "VH-AVH", "VH-nonAVH","AH-AVH","AH-nonAVH","AVH-nonAVH")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
covars_str <- c('raz', 'age', 'sex', 'T1_site')
str<-names(bsnip)[grep("Left_V1|Left_V2|Right_V1|Right_V2|Left_MT|Right_MT|Left_fusiform|Right_fusiform|Left_lingual|Right_lingual|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'avh', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'avh', covars_r)
results_str[results_str$p_adj<0.2,]


#compare HC to proband quantiles
bsnip$hcquant<- NA
bsnip$hcquant[bsnip$phenotype == "HC"] <- 'HC'
bsnip$hcquant[bsnip$phenotype == "Probands" & bsnip$quantile == "1"] <- 'q1'
bsnip$hcquant[bsnip$phenotype == "Probands" & bsnip$quantile == "2"] <- 'q2'
bsnip$hcquant[bsnip$phenotype == "Probands" & bsnip$quantile == "3"] <- 'q3'
bsnip$hcquant[bsnip$phenotype == "Probands" & bsnip$quantile == "4"] <- 'q4'
bsnip$hcquant <- as.factor(as.character(bsnip$hcquant))
table(bsnip$quantile)
table(bsnip$hcquant)

#run contrast for HC to proband quantiles in area 
bsnip$vh <- factor(bsnip$vh, levels = c('HC', 'q1', 'q2', 'q3', 'q4'))
pheno.mat<- rbind(c(1,-1,0,0,0), c(1,0,-1,0,0), c(1,0,0,-1,0),c(1,0,0,0,-1))
rownames(pheno.mat)<- c("HC-q1", "HC-q2", "HC-q3","HC-q4")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
covars_str <- c('raz', 'age', 'sex', 'T1_site')
str<-names(bsnip)[grep("Left_V1_exvivo_area|Left_V2_exvivo_area|Right_V1_exvivo_area|Right_V2_exvivo_area
|Left_MT_exvivo_area|Right_MT_exvivo_area|Left_fusiform|Right_fusiform_area|Left_lingual_area|Right_lingual_area|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'hcquant', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'hcquant', covars_r)
results_str[results_str$p_adj<0.2,]

#run contrast for HC to proband quantiles in volume
bsnip$vh <- factor(bsnip$vh, levels = c('HC', 'q1', 'q2', 'q3', 'q4'))
pheno.mat<- rbind(c(1,-1,0,0,0), c(1,0,-1,0,0), c(1,0,0,-1,0),c(1,0,0,0,-1))
rownames(pheno.mat)<- c("HC-q1", "HC-q2", "HC-q3","HC-q4")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
covars_str <- c('raz', 'age', 'sex', 'T1_site')
str<-names(bsnip)[grep("Left_V1_exvivo_volume|Left_V2_exvivo_volume|Right_V1_exvivo_volume|Right_V2_exvivo_volume
|Left_MT_exvivo_volume|Right_MT_exvivo_volume|Left_fusiform_exvivo_volume|Right_fusiform_volume|Left_lingual_volume|Right_lingual_volume|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'hcquant', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'hcquant', covars_r)
results_str[results_str$p_adj<0.2,]

#run contrast for HC to proband quantiles in thickness
bsnip$vh <- factor(bsnip$vh, levels = c('HC', 'q1', 'q2', 'q3', 'q4'))
pheno.mat<- rbind(c(1,-1,0,0,0), c(1,0,-1,0,0), c(1,0,0,-1,0),c(1,0,0,0,-1))
rownames(pheno.mat)<- c("HC-q1", "HC-q2", "HC-q3","HC-q4")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
covars_str <- c('raz', 'age', 'sex', 'T1_site')
str<-names(bsnip)[grep("Left_V1_exvivo_thickness|Left_V2_exvivo_thickness|Right_V1_exvivo_thickness|Right_V2_exvivo_thickness
|Left_MT_exvivo_thickness|Right_MT_exvivo_thickness|Left_fusiform_thickness|Right_fusiform_thickness|Left_lingual_thickness|Right_lingual_thickness|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'hcquant', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'hcquant', covars_r)
results_str[results_str$p_adj<0.1,]


#compare HC to proband quantiles
bsnip$proquant<- NA
bsnip$proquant[bsnip$phenotype == "Probands"] <- 'Proband'
bsnip$proquant[bsnip$phenotype == "HC" & bsnip$quantile == "1"] <- 'q1'
bsnip$proquant[bsnip$phenotype == "HC" & bsnip$quantile == "2"] <- 'q2'
bsnip$proquant[bsnip$phenotype == "HC" & bsnip$quantile == "3"] <- 'q3'
bsnip$proquant[bsnip$phenotype == "HC" & bsnip$quantile == "4"] <- 'q4'
bsnip$proquant <- as.factor(as.character(bsnip$proquant))
table(bsnip$phenotype, bsnip$quantile)
table(bsnip$quantile)
table(bsnip$proquant)
str(bsnip$proquant)

#run contrast for Proband to HC quantiles in area 
bsnip$vh <- factor(bsnip$vh, levels = c('Pro', 'q1', 'q2', 'q3', 'q4'))
pheno.mat<- rbind(c(1,-1,0,0,0), c(1,0,-1,0,0), c(1,0,0,-1,0),c(1,0,0,0,-1))
rownames(pheno.mat)<- c("Pro-q1", "Pro-q2", "Pro-q3","Pro-q4")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
covars_str <- c('raz', 'age', 'sex', 'T1_site')
str<-names(bsnip)[grep("Left_V1_exvivo_area|Left_V2_exvivo_area|Right_V1_exvivo_area|Right_V2_exvivo_area
|Left_MT_exvivo_area|Right_MT_exvivo_area|Left_fusiform|Right_fusiform_area|Left_lingual_area|Right_lingual_area|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'proquant', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'proquant', covars_r)
results_str[results_str$p_adj<0.1,]

#run contrast for Proband to HC quantiles in volume 
bsnip$vh <- factor(bsnip$vh, levels = c('Pro', 'q1', 'q2', 'q3', 'q4'))
pheno.mat<- rbind(c(1,-1,0,0,0), c(1,0,-1,0,0), c(1,0,0,-1,0),c(1,0,0,0,-1))
rownames(pheno.mat)<- c("Pro-q1", "Pro-q2", "Pro-q3","Pro-q4")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
covars_str <- c('raz', 'age', 'sex', 'T1_site')
str<-names(bsnip)[grep("Left_V1_exvivo_volume|Left_V2_exvivo_volume|Right_V1_exvivo_volume|Right_V2_exvivo_volume
|Left_MT_exvivo_volume|Right_MT_exvivo_volume|Left_fusiform_exvivo_volume|Right_fusiform_volume|Left_lingual_volume|Right_lingual_volume|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'proquant', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'proquant', covars_r)
results_str[results_str$p_adj<0.1,]

bsnip$vh <- factor(bsnip$vh, levels = c('HC', 'q1', 'q2', 'q3', 'q4'))
pheno.mat<- rbind(c(1,-1,0,0,0), c(1,0,-1,0,0), c(1,0,0,-1,0),c(1,0,0,0,-1))
rownames(pheno.mat)<- c("HC-q1", "HC-q2", "HC-q3","HC-q4")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
covars_str <- c('raz', 'age', 'sex', 'T1_site')
str<-names(bsnip)[grep("Left_V1_exvivo_thickness|Left_V2_exvivo_thickness|Right_V1_exvivo_thickness|Right_V2_exvivo_thickness
|Left_MT_exvivo_thickness|Right_MT_exvivo_thickness|Left_fusiform_thickness|Right_fusiform_thickness|Left_lingual_thickness|Right_lingual_thickness|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'hcquant', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'hcquant', covars_r)
results_str[results_str$p_adj<0.2,]
#run contrast for Proband to HC quantiles in thickness
bsnip$vh <- factor(bsnip$vh, levels = c('Pro', 'q1', 'q2', 'q3', 'q4'))
pheno.mat<- rbind(c(1,-1,0,0,0), c(1,0,-1,0,0), c(1,0,0,-1,0),c(1,0,0,0,-1))
rownames(pheno.mat)<- c("Pro-q1", "Pro-q2", "Pro-q3","Pro-q4")
covars_r<-c('raz','site', 'age', 'sex', 'eTIV')
covars_str <- c('raz', 'age', 'sex', 'T1_site','EstimatedTotalIntraCranialVol')
str<-names(bsnip)[grep("Left_V1_exvivo_thickness|Left_V2_exvivo_thickness|Right_V1_exvivo_thickness|Right_V2_exvivo_thickness
|Left_MT_exvivo_thickness|Right_MT_exvivo_thickness|Left_fusiform_thickness|Right_fusiform_thickness|Left_lingual_thickness|Right_lingual_thickness|
  Left.Thalamus|Right.Thalamus|Left.Hippocampus|Right.Hippocampus", names(bsnip))]
results_str<- RunMat(str, pheno.mat, bsnip, 'proquant', covars_r)
results_str$p_adj<-p.adjust(p=results_str$p, method = 'fdr')
results_str<-get_runcontrast_d(bsnip,results_str, 'proquant', covars_r)
results_str[results_str$p_adj<0.2,]

quantile <- c('quantile')
ctq<- c('ctq_total', 'ctqscore_ea', 'ctqscore_sa', 'ctqscore_pa', 'ctqscore_en', 'ctqscore_pn')
covars_r<-c('raz','age', 'sex', 'site','eTIV')
covars_str <- c('raz', 'age', 'sex', 'T1_site','EstimatedTotalIntraCranialVol')

#run correlation with HC and area
area<- c('Left_V1_exvivo_area', 'Left_V2_exvivo_area', 'Right_V1_exvivo_area', 'Right_V2_exvivo_area',
         'Left_MT_exvivo_area', 'Right_MT_exvivo_area', 'Left_fusiform_area', 'Right_fusiform_area', 'Left_lingual_area', 'Right_lingual_area',
         'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
vc_area<- get_correlations(df = bsnip, meas1 = ctq, meas2 = area, covars1 = covars_r, covars2 = covars_str, method = 'spearman')
vc_area[vc_area$p_value<0.1,]

#run correlation with HC and area
area<- c('Left_V1_exvivo_area', 'Left_V2_exvivo_area', 'Right_V1_exvivo_area', 'Right_V2_exvivo_area',
'Left_MT_exvivo_area', 'Right_MT_exvivo_area', 'Left_fusiform_area', 'Right_fusiform_area', 'Left_lingual_area', 'Right_lingual_area',
  'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
vc_area<- get_correlations(df = HC, meas1 = ctq, meas2 = area, covars1 = covars_r, covars2 = covars_str, method = 'spearman')
vc_area[vc_area$p_value<0.1,]

#run correlation with HC and volume
vol<- c('Left_V1_exvivo_volume', 'Left_V2_exvivo_volume', 'Right_V1_exvivo_volume', 'Right_V2_exvivo_volume',
        'Left_MT_exvivo_volume', 'Right_MT_exvivo_volume', 'Left_fusiform_volume', 'Right_fusiform_volume', 'Left_lingual_volume', 'Right_lingual_volume',
        'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
vc_vol<- get_correlations(HC, meas1 = ctq, meas2 = vol, covars1 = covars_r, covars2 = covars_str, method = 'spearman')
vc_vol[vc_vol$adjusted_p_value<1,]

#run correlation with HC and thickness
thick<- c('Left_V1_exvivo_thickness', 'Left_V2_exvivo_thickness', 'Right_V1_exvivo_thickness', 'Right_V2_exvivo_thickness',
        'Left_MT_exvivo_thickness', 'Right_MT_exvivo_thickness', 'Left_fusiform_thickness', 'Right_fusiform_thickness', 'Left_lingual_thickness', 'Right_lingual_thickness',
        'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
vc_thick<- get_correlations(HC, meas1 = ctq, meas2 = thick, covars1 = covars_str, covars2 = covars_r, method = 'spearman')
vc_thick[vc_thick$adjusted_p_value<1,]

#run correlation with nonVH and area
area<- c('Left_V1_exvivo_area', 'Left_V2_exvivo_area', 'Right_V1_exvivo_area', 'Right_V2_exvivo_area',
         'Left_MT_exvivo_area', 'Right_MT_exvivo_area', 'Left_fusiform_area', 'Right_fusiform_area', 'Left_lingual_area', 'Right_lingual_area',
         'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
vc_area<- get_correlations(nonVH, meas1 = ctq, meas2 = area, covars1 = covars_str, covars2 = covars_r, method = 'spearman')
vc_area[vc_area$adjusted_p_value<1,]

#run correlation with nonVH and volume
vol<- c('Left_V1_exvivo_volume', 'Left_V2_exvivo_volume', 'Right_V1_exvivo_volume', 'Right_V2_exvivo_volume',
        'Left_MT_exvivo_volume', 'Right_MT_exvivo_volume', 'Left_fusiform_volume', 'Right_fusiform_volume', 'Left_lingual_volume', 'Right_lingual_volume',
        'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
#thalamus and hippocampus are volume
vc_vol<- get_correlations(nonVH, meas1 = ctq, meas2 = vol, covars1 = covars_str, covars2 = covars_r, method = 'spearman')
vc_vol[vc_vol$adjusted_p_value<1,]

#run correlation with nonVH and thickness
thick<- c('Left_V1_exvivo_thickness', 'Left_V2_exvivo_thickness', 'Right_V1_exvivo_thickness', 'Right_V2_exvivo_thickness',
          'Left_MT_exvivo_thickness', 'Right_MT_exvivo_thickness', 'Left_fusiform_thickness', 'Right_fusiform_thickness', 'Left_lingual_thickness', 'Right_lingual_thickness',
          'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
vc_thick<- get_correlations(nonVH, meas1 = ctq, meas2 = thick, covars1 = covars_str, covars2 = covars_r, method = 'spearman')
vc_thick[vc_thick$adjusted_p_value<1,]

#run correlation with VH and area
area<- c('Left_V1_exvivo_area', 'Left_V2_exvivo_area', 'Right_V1_exvivo_area', 'Right_V2_exvivo_area',
         'Left_MT_exvivo_area', 'Right_MT_exvivo_area', 'Left_fusiform_area', 'Right_fusiform_area', 'Left_lingual_area', 'Right_lingual_area',
         'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
vc_area<- get_correlations(VH, meas1 = ctq, meas2 = area, covars1 = covars_str, covars2 = covars_r, method = 'spearman')
vc_area[vc_area$adjusted_p_value<0.1,]

#run correlation with VH and volume
vol<- c('Left_V1_exvivo_volume', 'Left_V2_exvivo_volume', 'Right_V1_exvivo_volume', 'Right_V2_exvivo_volume',
        'Left_MT_exvivo_volume', 'Right_MT_exvivo_volume', 'Left_fusiform_volume', 'Right_fusiform_volume', 'Left_lingual_volume', 'Right_lingual_volume',
        'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
vc_vol<- get_correlations(VH, meas1 = ctq, meas2 = vol, covars1 = covars_str, covars2 = covars_r, method = 'spearman')
vc_vol[vc_vol$adjusted_p_value<0.1,]

#run correlation with VH and thickness
thick<- c('Left_V1_exvivo_thickness', 'Left_V2_exvivo_thickness', 'Right_V1_exvivo_thickness', 'Right_V2_exvivo_thickness',
          'Left_MT_exvivo_thickness', 'Right_MT_exvivo_thickness', 'Left_fusiform_thickness', 'Right_fusiform_thickness', 'Left_lingual_thickness', 'Right_lingual_thickness',
          'Left.Thalamus', 'Right.Thalamus', 'Left.Hippocampus', 'Right.Hippocampus')
vc_thick<- get_correlations(VH, meas1 = ctq, meas2 = thick, covars1 = covars_str, covars2 = covars_r, method = 'spearman')
vc_thick[vc_thick$adjusted_p_value<0.1,]

summary(aov(bsnip2$Left_MT_exvivo_thickness~bsnip2$quantile*bsnip2$phenotype+bsnip2$age+bsnip2$sex+bsnip2$raz+bsnip2$site))
x <- aov(bsnip2$quantile~bsnip2$Left_MT_exvivo_thickness+bsnip2$age+bsnip2$sex+bsnip2$raz+bsnip2$site)
TukeyHSD(x)

#########structure
#mediation analysis
plot(bsnip$quantile, bsnip$Left.Hippocampus)
install.packages("gvlma")
install.packages("stargazer")
install.packages("bda")
install.packages("multilevel")
install.packages("MBESS")
library(MBESS)
library(mediation)
library(rockchalk) #Graphing simple slopes; moderation
library(multilevel) #Sobel Test
library(bda) #Another Sobel Test option
library(gvlma) #Testing Model Assumptions
library(stargazer) #Handy regression tables
#X = CTQ, M = structure, Y = VH
class(bsnip$ctq_total)
class(bsnip$vh)
bsnip$vh = as.numeric(bsnip$vh)
class(bsnip$vh)
proband$vh2= as.numeric(proband$vh2)
#test if Y and X have a relationship
model.0 <- lm(vh~ctq_total, data = bsnip)
summary(model.0)
#test if M and X have a relationship - average right and left
model.M <- lm(Left_MT_exvivo_thickness~ctq_total, bsnip)
summary(model.M)
#test effect of X on Y through Mmodel.Y <- lm(vh~Left_MT_exvivo_thickness+ctq_total, bsnip)
summary(model.Y)
#4. Reversed Path C (Y on X, controlling for M)
fitc <- lm(ctq_total ~ vh + Left_MT_exvivo_thickness, bsnip)
summary(fitc)
#summary table
stargazer(model.0, model.M, model.Y, fitc, type = "text", title = "Baron and Kenny Method")
#interpreting Baron and Kenny
mediation.test(bsnip$Left_V1_exvivo_area,bsnip$ctq_total,bsnip$vh)
?sobel
#sobel(Meddata$X, Meddata$M, Meddata$Y)
sobel(bsnip$ctq_total,bsnip$Left_V1_exvivo_area,bsnip$vh)
proband$vh2 = as.numeric(proband$vh2)
class(proband$vh2)

thick_adj<- c('Left_V1_exvivo_thickness', 'Left_V2_exvivo_thickness', 'Right_V1_exvivo_thickness', 'Right_V2_exvivo_thickness',
              'Left_MT_exvivo_thickness', 'Right_MT_exvivo_thickness')
area_adj<- c('Left_V1_exvivo_area', 'Left_V2_exvivo_area', 'Right_V1_exvivo_area', 'Right_V2_exvivo_area',
         'Left_MT_exvivo_area', 'Right_MT_exvivo_area')
bsnip[, thick_adj]=sapply(thick_adj, function(x) returnAdj(data = bsnip, measure = x, covars = c('age','sex','raz','site')))
bsnip[, area_adj]=sapply(area_adj, function(x) returnAdj(data = bsnip, measure = x, covars = c('age','sex','raz','site', 'EstimatedTotalIntraCranialVol')))

#MBESS
?mediation
bsnip$vh<- NA
bsnip$vh[bsnip$p3_anhal_b == "1"] = "VH"
bsnip$vh[bsnip$group == "NC"] = "HC"
bsnip$vh[bsnip$p3_anhal_b == "0"&bsnip$group!="NC"] = "nonVH"
class(bsnip$vh)
bsnip$vh2[bsnip$p3_anhal_b == "1"] = "VH"
bsnip$vh2[bsnip$p3_anhal_b == "0"&bsnip$group!="NC"] = "nonVH"
#X = CTQ, M = structure, Y = VH
class(bsnip$ctq_total)
#bsnip$vh[bsnip$p3_anhal_b == "1"] = "VH"
#bsnip$vh[bsnip$p3_anhal_b == "0"&bsnip$group!="NC"] = "nonVH"
#bsnip$vh[bsnip$group == "NC"] = "HC"
bsnip$vh[bsnip$p3_anhal_b == "1"] = "2"
bsnip$vh[bsnip$p3_anhal_b == "0"&bsnip$group!="NC"] = "1"
bsnip$vh[bsnip$group == "NC"] = "0"
bsnip$vh2[bsnip$p3_anhal_b == "1"] = "1"
bsnip$vh2[bsnip$p3_anhal_b == "0"&bsnip$group!="NC"] = "0"
bsnip$vh
bsnip$vh <- as.numeric(bsnip$vh)
bsnip$vh2 <- as.numeric(bsnip$vh2)
class(bsnip$vh2)
#creating average structure variable
bsnip$avgmt<- (bsnip$Left_MT_exvivo_thickness + bsnip$Right_MT_exvivo_thickness)
bsnip$avgmt <- bsnip$avgmt/2
bsnip$avgv1<- (bsnip$Left_V1_exvivo_thickness+ bsnip$Right_V1_exvivo_thickness)
bsnip$avgv2<- (bsnip$Left_V2_exvivo_thickness+ bsnip$Right_V2_exvivo_thickness)
bsnip$avgv2<- bsnip$avgv2/2
bsnip$avgv2area<- (bsnip$Left_V2_exvivo_area+ bsnip$Right_V2_exvivo_area)
bsnip$avgv2area<- bsnip$avgv2area/2
bsnip$avgv1area<- (bsnip$Left_V1_exvivo_area+ bsnip$Right_V1_exvivo_area)
bsnip$avgv1area<- bsnip$avgv1area/2
bsnip$avgstr
#testing indirect effect
mediation(bsnip$ctq_total,bsnip$avgstr,bsnip$vh,bootstrap=TRUE,
which.boot="both",B=1000,conf.level=.95)
bsnip= bsnip[!is.na(bsnip$avgstr),]
#Mediate package
fit1<- lm(avgmt~ctq_total, data=bsnip)
fit2<- lm(vh~avgmt+ctq_total, data=bsnip)
fit1<- lm(avgv1~ctq_total, data=bsnip)
fit2<- lm(vh~avgv1+ctq_total, data=bsnip)
fit1<- lm(avgv2~ctq_total, data=bsnip)
fit2<- lm(vh~avgv2+ctq_total, data=bsnip)
fit1<- lm(avgv2area~ctq_total, data=bsnip)
fit2<- lm(vh~avgv2area+ctq_total, data=bsnip)
fit1<- lm(avgv1area~ctq_total, data=bsnip)
fit2<- lm(vh~avgv1area+ctq_total, data=bsnip)
fitMed <- mediate(fit1, fit2, treat="ctq_total", mediator="avgmt")
fitMed <- mediate(fit1, fit2, treat="ctq_total", mediator="avgv1")
fitMed <- mediate(fit1, fit2, treat="ctq_total", mediator="avgv2")
fitMed <- mediate(fit1, fit2, treat="ctq_total", mediator="avgv2area")
fitMed <- mediate(fit1, fit2, treat="ctq_total", mediator="avgv1area")

summary(fit1)
summary(fit2)
summary(fitMed)
gvlma(fit1)
gvlma(fit2)
#Bootstrap alternative
fitMedBoot <- mediate(fit1, fit2, boot=TRUE, sims=1000, treat="ctq_total", mediator="avgmt")
fitMedBoot <- mediate(fit1, fit2, boot=TRUE, sims=1000, treat="ctq_total", mediator="avgv1")
fitMedBoot <- mediate(fit1, fit2, boot=TRUE, sims=1000, treat="ctq_total", mediator="avgv2")
fitMedBoot <- mediate(fit1, fit2, boot=TRUE, sims=1000, treat="ctq_total", mediator="avgv2area")
fitMedBoot <- mediate(fit1, fit2, boot=TRUE, sims=1000, treat="ctq_total", mediator="avgv1area")
summary(fitMedBoot)
plot(fitMedBoot)




  
