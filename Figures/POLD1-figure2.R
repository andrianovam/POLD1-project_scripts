library(ggplot2)
library(ggtext)
library(dplyr)
library("MutationalPatterns")
library("this.path")
library(cowplot)
library(gridExtra)
library(reshape2)
library(ggpubr)

setwd(dirname(this.path()))

#Figure 2A

df = read.table("Input_files/Figure2A_fibros_mutrate_data.txt", header=T)

samples_levels = c("IV.5","IV.3", "IV.2", "III.2", "III.4", "IV.4")
df_SNV = df[,c(1,9,13)]
colnames(df_SNV)[2]<-"Number"
df_SNV$Type = "SNVs"
df_indels = df[,c(1,10,13)]
colnames(df_indels)[2]<-"Number"
df_indels$Type = "Indels"
df_indels$Number = df_indels$Number*5
df_all <- rbind(df_SNV, df_indels)

tiff(filename="Output_plots/Figure2A_SNVs_indels_fibros_merged.tiff", width=10, height=6, res=300, units='cm')

(ggplot(df_all) +
    geom_bar(aes(x = factor(Sample,levels = samples_levels), y = Number, alpha=Type, fill=Status), stat = "identity",position = "dodge") +
    scale_y_continuous(name = "SNVs",sec.axis = sec_axis( trans=~./5, name="Indels")) +
    theme_bw() +
    scale_fill_manual(values=c('orange','cornflowerblue')) +
    scale_alpha_manual(values=c(0.6,1)) +
    xlab("")+
    theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=7))
  
)
dev.off()
  
#Figure 2B (COSMIC Signatures in p40 - refitted subset of signatures by SigFit)

df = read.table("Input_files/Figure2B_fibros_mutrate_data_refitted.txt", header=T)
df_subset = df[,c(1,12,3,4,5,6,7,8)]
result_df = melt(df_subset)

sigs_levels = c("SBS10c","SBS36", "SBS37", "SBS45", "SBS93","SBS5")
samples_levels = c('IV.4', 'III.4', 'III.2', 'IV.2', 'IV.3', 'IV.5')

tiff(filename="Output_plots/Figure2B_signatures_fibros_p40_refitted.tiff", width=10, height=6, res=300, units='cm')

(ggplot(aes(y=value, x=factor(Sample,levels = samples_levels), fill = factor(variable, levels = sigs_levels)), data = result_df)
  + geom_bar(stat = "identity")
  + scale_y_continuous(labels = scales::percent) 
  + facet_grid(cols = vars(result_df$Status), scale = 'free_x',space = "free_x")
  + theme_bw() 
  + scale_fill_manual(values = c("darkred","bisque1","bisque3","darkgrey","beige","darkblue"),labels = c("SBS10c\n(POLD1 deficiency)", "SBS36", "SBS37", "SBS45\n(possible artefact)", "SBS93", "SBS5\n(clock-like)"))
  + xlab("")
  + ylab("")
  + theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=7))
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  + theme(strip.text = element_text(size = 10))
  + theme(legend.key.height = unit(0.5, "cm"),legend.margin=margin(t = 0, unit='cm'),legend.title=element_blank())
)

dev.off()

#Figure 2C-E - PCA for somatic mutations accumulated in the experiment

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

vcf_files <- list.files(path = './Input_files/Fibros_p40/', pattern="*.vcf", full.names=T)
vcf_files
sample_names <- c("III.2","III.4","IV.2","IV.3","IV.4", "IV.5")
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,  type = "all")
summary(vcfs)

#function that takes mutation matrix with counts created by mut_matrix function from MutationalPatterns and
#calculates 3-nucleotide context proportions within mutation type
create_mut_mat_freq_by_type <- function(mut_mat){
  mut_mat_CA = mut_mat[1:16,]
  mut_mat_CG = mut_mat[17:32,]
  mut_mat_CT = mut_mat[33:48,]
  mut_mat_TA = mut_mat[49:64,]
  mut_mat_TC = mut_mat[65:80,]
  mut_mat_TG = mut_mat[81:96,]
  mut_mat_CA_freq = apply(mut_mat_CA, 1, function(x) x/colSums(mut_mat_CA))
  mut_mat_CG_freq = apply(mut_mat_CG, 1, function(x) x/colSums(mut_mat_CG))
  mut_mat_CT_freq = apply(mut_mat_CT, 1, function(x) x/colSums(mut_mat_CT))
  mut_mat_TA_freq = apply(mut_mat_TA, 1, function(x) x/colSums(mut_mat_TA))
  mut_mat_TC_freq = apply(mut_mat_TC, 1, function(x) x/colSums(mut_mat_TC))
  mut_mat_TG_freq = apply(mut_mat_TG, 1, function(x) x/colSums(mut_mat_TG))
  mut_mat_freq_by_type = cbind(mut_mat_CA_freq,mut_mat_CG_freq,mut_mat_CT_freq,mut_mat_TA_freq,mut_mat_TG_freq,mut_mat_TC_freq)
  return(mut_mat_freq_by_type)
}

mut_mat_somatic <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
mut_mat_freq_by_type_somatic = create_mut_mat_freq_by_type(mut_mat_somatic)

Status<-c("L474P","L474P","L474P","wt","L474P","wt")
pol<-c("POLD","POLD","POLD", "Non-carier","POLD","Non-carier")

res.pca_by_type <-prcomp(mut_mat_freq_by_type_somatic)
result_pc_by_type<-as.data.frame(res.pca_by_type$x)
result_pc_by_type$PolD_status = Status
result_pc_by_type$Pol = pol

nudges_y = c(-0.01, 0, 0, 0, 0.03, 0.01)
nudges_x = c(0.025, 0.025, 0.025, 0.025, 0.01, 0.025)

#2C
tiff(filename="Output_plots/Figure2C_fibros_p40_pca.tiff", width=8, height=8, res=300, units='cm')
(ggplot(result_pc_by_type, aes(x=-PC1, y=PC2, col=PolD_status))
  + geom_point(size=2)
  + theme_bw()
  + geom_text(label=rownames(result_pc_by_type),nudge_x = nudges_x, nudge_y = nudges_y, check_overlap = T, size=3)
  + theme(legend.position = "bottom")+scale_color_manual(values=c('orange','cornflowerblue'), name="*POLD1* status")
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8))
  + theme(legend.margin =margin(r=0,l=0,t=-10,b=0))
  + theme(legend.box.margin=margin(0,0,-5,0))
  + xlim(-0.18,0.15)
  + ylim(-0.18,0.13)
  + theme(legend.title = element_markdown()))
dev.off()

#2E
weights_pc_by_type<-as.data.frame(res.pca_by_type$rotation)
weights_pc_by_type$Mut_type = rownames(weights_pc_by_type)
weights_pc_by_type$Mutation = substr(weights_pc_by_type$Mut_type,3,5)
weights_pc_by_type$minus_PC1 = -weights_pc_by_type$PC1

muttype_level = c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T" )
x_labels = c("A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T" )


tiff(filename="Output_plots/Figure2E_PC1_spectrum.tiff", width=18, height=6, res=300, units='cm')
(ggplot(weights_pc_by_type)
  + geom_bar(aes(x=factor(Mut_type,levels = muttype_level), y=minus_PC1, fill=Mutation),stat="identity")
  + theme_bw()
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "top",legend.box = "horizontal",strip.text = element_text(size = 10))
  + guides(fill = guide_legend(nrow = 1))
  + scale_x_discrete(labels=x_labels)
  + scale_fill_manual(values = c("deepskyblue2","black","red","grey","darkolivegreen3","lightpink"))
  + ylab("-PC1")
  + theme(legend.key.size = unit(0.4, "cm"))
  + theme(axis.text.y=element_text(size=10),axis.text.x=element_text(size=6),axis.title=element_text(size=10),legend.text=element_text(size=8))
  + xlab(""))
dev.off()

#2D
Tissue<-c("Fibros","Fibros","Fibros","Fibros_wt","Fibros","Fibros_wt")
result_pc_by_type$Tissue = Tissue
result_pc_by_type$Data = "Experiment"
#download data for crypts and project to pc space
vcf_files_colon_PolD <- list.files(path = './Input_files/crypts/', pattern="*.vcf", full.names = T)
vcfs_colon_PolD <- read_vcfs_as_granges(vcf_files_colon_PolD, vcf_files_colon_PolD, ref_genome,  type = "snv")
mut_mat_colon_PolD <- mut_matrix(vcf_list = vcfs_colon_PolD, ref_genome = ref_genome)
mut_mat_freq_by_type_colon_PolD = create_mut_mat_freq_by_type(mut_mat_colon_PolD)
colon_PolD_pca_by_type<-as.data.frame(predict(res.pca_by_type, mut_mat_freq_by_type_colon_PolD))
colon_PolD_pca_by_type$Tissue = "colon"
colon_PolD_pca_by_type$Pol = "PolD"
colon_PolD_pca_by_type$PolD_status = c(rep("S478N", 35), rep("L474P", 8), rep("D316N",9))
colon_PolD_pca_by_type$Data = "Accumulated with age"
#download data for sperm and project to pc space
vcf_files_duplex_PolD <- list.files(path = './Input_files/blood_sperm/', pattern="*.vcf", full.names = T)
vcfs_duplex_PolD <- read_vcfs_as_granges(vcf_files_duplex_PolD, vcf_files_duplex_PolD, ref_genome,  type = "snv")
mut_mat_duplex_PolD <- mut_matrix(vcf_list = vcfs_duplex_PolD, ref_genome = ref_genome)
mut_mat_freq_by_type_duplex_PolD = create_mut_mat_freq_by_type(mut_mat_duplex_PolD)
duplex_PolD_pca_by_type<-as.data.frame(predict(res.pca_by_type, mut_mat_freq_by_type_duplex_PolD))
duplex_PolD_pca_by_type$Tissue = c("sperm","blood")
duplex_PolD_pca_by_type$Pol = "PolD"
duplex_PolD_pca_by_type$PolD_status = "S478N"
duplex_PolD_pca_by_type$Data = "Accumulated with age"

#download data for fibros p0
vcf_files_fibros_p0 <- list.files(path = './Input_files/Fibros_p0/', pattern="*.vcf", full.names = T)
vcfs_fibros_p0 <- read_vcfs_as_granges(vcf_files_fibros_p0, vcf_files_fibros_p0, ref_genome,  type = "snv")
mut_mat_fibros_p0 <- mut_matrix(vcf_list = vcfs_fibros_p0, ref_genome = ref_genome)
mut_mat_freq_by_type_fibros_p0 = create_mut_mat_freq_by_type(mut_mat_fibros_p0)
fibros_p0_pca_by_type<-as.data.frame(predict(res.pca_by_type, mut_mat_freq_by_type_fibros_p0))
fibros_p0_pca_by_type$Tissue = c("Fibros_p0","Fibros_p0","Fibros_p0","Fibros_p0","Fibros_p0_wt","Fibros_p0","Fibros_p0_wt","Fibros_p0")
fibros_p0_pca_by_type$Pol = "PolD"
fibros_p0_pca_by_type$PolD_status = c("L474P","L474P","L474P","L474P","wt","L474P","wt","L474P")
fibros_p0_pca_by_type$Data = "Accumulated with age"

#merge results for fibros, crypts and sperm in one dataframe
result_by_type = rbind(result_pc_by_type, colon_PolD_pca_by_type)
result_by_type = rbind(result_by_type, duplex_PolD_pca_by_type)
result_by_type = rbind(result_by_type, fibros_p0_pca_by_type)

x_levels = c("Fibro_wt", "Fibro_PolD", "colon","blood", "sperm")
tissue_levels = c("Fibros_wt", "Fibros", "Fibros_p0_wt", "Fibros_p0", "colon", "blood","sperm")
mutation_levels = c("wt", "L474P","D316N", "S478N")

tiff(filename="Output_plots/Figure2D_other_tissues_PC1.tiff", width=9, height=6, res=300, units='cm')
(ggplot(result_by_type, aes(x=factor(Tissue, levels = tissue_levels), y= -PC1))
  + geom_boxplot()
  + geom_point(aes(col=PolD_status), size=1)
  + theme_bw()
  + xlab("")
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8))
  + theme(legend.margin =margin(r=0,l=20,t=-35,b=0))
  + theme(legend.box.margin=margin(0,20,-35,0))
  + theme(legend.spacing.x = unit(0.05, 'cm'))
  + theme(legend.position = "bottom",legend.box="vertical",legend.margin=margin(0,0,30,0))
  + ylim(-0.2,0.35)
  + scale_x_discrete(labels=c("Fibros_p0_wt" = "Fibros\nwt\n(2)", "Fibros_p0" = "Fibros\nL474P\n(8)" ,"Fibros_wt" = "Fibros\nwt\n(2)", "Fibros" = "Fibros\nL474P\n(4)","colon" = "colon\n(52)", "blood"="blood\n(1)","sperm"="sperm\n(1)"))
  + scale_color_discrete(name = "Mutation")
  + facet_grid(cols=vars(factor(result_by_type$Data, levels = c("Experiment","Accumulated with age"))),scales="free_x", space="free_x"))
dev.off()


#Figure 2F number of mutations germline experiment vs public trios

df = read.table("Input_files/public_trio_data.txt")
head(df)
colnames(df) = c("chr","coord","mut","family_id","age","dataset")
count_mut_per_family = aggregate(mut ~ family_id, data = df[,c(3,4)], FUN = length)
head(count_mut_per_family)
count_mut_per_family$POLD = "Published trios"
our_data = data.frame(family_id = c("IV.1","IV.2","IV.3","IV.4","IV.5","IV.6","IV.7","IV.8"), mut=c(64,62,71,61,64,88,24,50), POLD = c("Mother","Mother","Mother","Father","Father","Father","Father","Both wt"))
data = rbind(count_mut_per_family, our_data)

levels_pold = c("Father","Mother", "Both wt", "Published trios")
levels_pold_our = c("Father","Mother", "Both wt")

tiff(filename="Output_plots/Figure2F_germline_counts.tiff", width=9, height=9, res=300, units='cm')

(ggplot()
  + geom_violin(data = data, aes(x=factor(POLD,levels = levels_pold),y=mut, fill=POLD), alpha=0.6)
  + geom_jitter(data = count_mut_per_family, aes(x=factor(POLD,levels = levels_pold),y=mut), size=0.05, alpha=0.3)
  + geom_violin(data = data, aes(x=factor(POLD,levels = levels_pold),y=mut, fill=POLD), alpha=0.6)
  + geom_point(data = our_data, aes(x=factor(POLD,levels = levels_pold_our),y=mut), size=0.5)
  + theme_bw()
  + theme(axis.title.x = element_markdown())
  + theme(axis.title.y = element_markdown())
  + ylab("Number of *de novo* mutations")+xlab("Parents *POLD1* status")
  + scale_fill_manual(values=c("indianred3","khaki3","mediumseagreen"))
  + scale_color_manual(values=c('black', "indianred4","khaki4","darkgreen"))
  + scale_x_discrete(labels=c("Father" = "Father\nL474P\n(4)", "Mother" = "Mother\nL474P\n(3)", "Both wt" = "Both wt\nour data\n(1)", "Published trios" = "Both wt\npublished\n(6241)"))
  + stat_summary(data = data, aes(x=factor(POLD,levels = levels_pold),y=mut), fun.data=mean_sdl, geom="pointrange",size=0.1, col='black')
  + theme(legend.position='none')
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8))
)

dev.off()

#Figure 2G percentage of POLD1 contexts germline experiment vs public trios

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

normal_trios = read.table("Input_files/public_trio_data_matrix.txt", header=T)
vcf_files_germline <- list.files(path = "Input_files/de_novo/", full.names = T)
vcf_files_germline
sample_names_germline <- c("IV.1", "IV.2","IV.3","IV.4", "IV.5","IV.6","IV.7", "IV.8")
status_germline<-c("Mother","Mother","Mother","Father","Father","Father","Father","Both wt", rep("Published trios", ncol(normal_trios)), "Father")
set_germline<-c("Main","Main","Main","Main","Main","Main","Main","Main", rep("Main", ncol(normal_trios)), "Add")
vcfs_germline <- read_vcfs_as_granges(vcf_files_germline, sample_names_germline, ref_genome,  type = "all")
summary(vcfs_germline)

mut_mat_germline <- mut_matrix(vcf_list = vcfs_germline, ref_genome = ref_genome)

sperm = read.table("Input_files/duplex_trinucl.csv", sep=";", header=T, row.names=1)

mut_mat_all = cbind(mut_mat_germline, normal_trios)
mut_mat_all = cbind(mut_mat_all, sperm$POLD1_sperm)
four_contexts = c("T[C>A]T", 'C[C>A]T', 'A[T>A]T','C[T>G]T')
subset_4contexts = mut_mat_all[unlist(four_contexts), ]
result = data.frame(sample = colnames(mut_mat_all), All_snv = colSums(mut_mat_all), four_context_prop = colSums(subset_4contexts), POLD = status_germline, Set = set_germline)
result$four_contexts_proportion = result$four_context_prop/result$All_snv * 100


levels_pold = c("Father","Mother", "Both wt", "Published trios")
levels_pold_our = c("Father","Mother", "Both wt")

data1 <- result %>% filter(Set=='Main')
data2 <- result %>% filter(Set =='Add')
our <- result %>% filter(POLD!='Published trios'&Set=='Main')
public <- result %>% filter(POLD =='Published trios'&Set=='Main')

tiff(filename="Output_plots/Figure2G_germline_4contexts_prop.tiff", width=8, height=8, res=300, units='cm')
(ggplot()
  + geom_violin(data = data1, aes(x=factor(POLD,levels = levels_pold),y=four_contexts_proportion, fill=POLD), alpha=0.6)
  + geom_jitter(data = public, aes(x=factor(POLD,levels = levels_pold),y=four_contexts_proportion), size=0.05, alpha=0.3)
  + geom_violin(data = data1, aes(x=factor(POLD,levels = levels_pold),y=four_contexts_proportion, fill=POLD), alpha=0.6)
  + geom_point(data = our, aes(x=factor(POLD,levels = levels_pold_our),y=four_contexts_proportion), size=0.5)
  + geom_point(data = data2, aes(x=factor(POLD,levels = levels_pold),y=four_contexts_proportion), size=1, col="indianred4")
  + geom_text(data=data2,aes(x=factor(POLD,levels = levels_pold),y=four_contexts_proportion, label="Sperm"),hjust=-0.2, size=3)
  + theme_bw()
  + theme(axis.title.x = element_markdown())
  + theme(axis.title.y = element_markdown())
  + ylab("Percent of *POLD1* contexts")+xlab("Parents *POLD1* status")
  + scale_fill_manual(values=c("indianred3","khaki3","mediumseagreen"))
  + scale_x_discrete(labels=c("Father" = "Father\nL474P\n(4)", "Mother" = "Mother\nL474P\n(3)", "Both wt" = "Both wt\nour data\n(1)", "Published trios" = "Both wt\npublished\n(6241)"))
  + stat_summary(data = data1, aes(x=factor(POLD,levels = levels_pold),y=four_contexts_proportion), fun.data=mean_sdl, geom="pointrange", color="black",size=0.2)
  + theme(legend.position='none')
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8))
)

dev.off()

#Figure 2H PC1 germline experiment vs public trios

normal_trios_by_type = read.table("Input_files/public_germline_pca_by_type.txt", header=T)
normal_trios_by_type$Status = "Published trios"

our_trios_by_type = read.table("Input_files/our_germline_pca_by_type.txt", header=T)
our_trios_by_type$Status = c("Father", "Father","Father","Mother","Father","Mother","Mother", "Both wt")

sperm_by_type = read.table("Input_files/sperm_duplex_germline_pca_freq_by_type.txt", header = T)
sperm_by_type = sperm_by_type[,1:6]
sperm_by_type$Status = "Father"
sperm_by_type$Set = "Add"

result_by_type = rbind(normal_trios_by_type, our_trios_by_type)
result_by_type$Set = "Main"
result_by_type = rbind(result_by_type, sperm_by_type[4,])

levels_pold = c("Father","Mother", "Both wt", "Published trios")
levels_pold_our = c("Father","Mother", "Both wt")

data1 <- result_by_type %>% filter(Set=='Main')
data2 <- result_by_type %>% filter(Set =='Add')
our <- result_by_type %>% filter(Status!='Published trios'&Set=='Main')
public <- result_by_type %>% filter(Status =='Published trios'&Set=='Main')

tiff(filename="Output_plots/Figure2H_germline_PC1_by_type.tiff", width=8, height=8, res=300, units='cm')
(ggplot()
  + geom_violin(data = data1, aes(x=factor(Status,levels = levels_pold),y=-PC1, fill=Status), alpha=0.6)
  + geom_jitter(data = public, aes(x=factor(Status,levels = levels_pold),y=-PC1), size=0.05, alpha=0.3)
  + geom_violin(data = data1, aes(x=factor(Status,levels = levels_pold),y=-PC1, fill=Status), alpha=0.6)
  + geom_point(data = our, aes(x=factor(Status,levels = levels_pold_our),y=-PC1), size=0.5)
  + geom_point(data = data2, aes(x=factor(Status,levels = levels_pold),y=-PC1), size=1, col="indianred4")
  + geom_text(data=data2,aes(x=factor(Status,levels = levels_pold),y=-PC1, label="Sperm"),hjust=-0.2, size=3)
  + theme_bw()
  + theme(axis.title.x = element_markdown())
  + theme(axis.title.y = element_markdown())
  + ylab("-PC1")+xlab("Parents *POLD1* status")
  + scale_fill_manual(values=c("indianred3","khaki3","mediumseagreen"))
  + scale_x_discrete(labels=c("Father" = "Father\nL474P\n(4)", "Mother" = "Mother\nL474P\n(3)", "Both wt" = "Both wt\nour data\n(1)", "Published trios" = "Both wt\npublished\n(6241)"))
  + stat_summary(data = data1, aes(x=factor(Status,levels = levels_pold),y=-PC1), fun.data=mean_sdl, geom="pointrange", color="black",size=0.2)
  + theme(legend.position='none')
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8))
)

dev.off()


