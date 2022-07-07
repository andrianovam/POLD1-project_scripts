library(ggplot2)
library(ggtext)
library(dplyr)
library("MutationalPatterns")
library("this.path")
library(cowplot)
library(gridExtra)
library(reshape2)

setwd(dirname(this.path()))

#Figure 2a

df = read.table("Input_files/Figure2A_signatures_fibros_p0.txt", header=T)
result_df = melt(df)

sigs_levels = c("SBS58","SBS10c","SBS1", "SBS5", "SBS7a", "SBS7b", "SBS7d")
samples_levels = c('III.2', 'III.4','IV.1','IV.2', 'IV.4','IV.6', 'IV.3', 'IV.5')

tiff(filename="Output_plots/Figure2a_signatures_fibros_p0.tiff", width=10, height=6, res=300, units='cm')

(ggplot(aes(y=value, x=factor(Samples,levels = samples_levels), fill = factor(variable, levels = sigs_levels)), data = result_df)
  + geom_bar(stat = "identity")
  + facet_grid(cols = vars(result_df$Status), scale = 'free_x',space = "free_x")
  + theme_bw() 
  + scale_fill_manual(values = c("grey", "darkred", "blue", "darkblue","darkseagreen1", "darkseagreen2", "darkseagreen3"), labels = c("SBS58\n(possible artefacts)", "SBS10c\n(POLD1 deficiency)", "SBS1\n(clock-like)", "SBS5\n(clock-like)", "SBS7a\n(UV radiation)", "SBS7b\n(UV radiation)", "SBS7d\n(UV radiation)"))
  + xlab("")
  + ylab("")
  + theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=7))
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  + theme(strip.text = element_text(size = 10))
  + theme(legend.key.height = unit(0.5, "cm"),legend.margin=margin(t = 0, unit='cm'),legend.title=element_blank())
)
   

dev.off()

#Figure 2b

df = read.table("Input_files/Figure2BC_fibros_mutrate_data.txt", header=T)

samples_levels = c("IV.5","IV.3", "IV.2", "III.2", "III.4", "IV.4")

tiff(filename="Output_plots/Figure2b_mutrate_fibros.tiff", width=7.5, height=7, res=300, units='cm')
(ggplot(df, aes(x=factor(Sample,levels = samples_levels), y=N_mut_all, fill=Status))
  + geom_bar(stat = 'identity',col="black")
  + theme_bw()
  + xlab("")
  + ylab("SNVs")
  + scale_fill_manual(values=c("darkgrey","white"))
  + coord_flip()
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8))
  + theme(legend.key.size = unit(0.4, "cm"), legend.position = "top")
  + theme(plot.margin = margin(0,0.5,0,0, "cm")))

dev.off()


#Figure 2c

df = read.table("Input_files/Figure2BC_fibros_mutrate_data.txt", header=T)

df_subset = df[,c(1,3,4,5,6,7,8)]
df_subset_prop = as.data.frame(apply(df_subset[,c(2:7)], 2, function(x) x/rowSums(df_subset[,c(2:7)])))
df_subset_prop$Sample = df_subset$Sample
df_subset_prop$Type = df$Status
result_df = melt(df_subset_prop)

sigs_levels = c("SBS10c","SBS36", "SBS37", "SBS45", "SBS93","SBS5")
samples_levels = c('IV.4', 'III.4', 'III.2', 'IV.2', 'IV.3', 'IV.5')

tiff(filename="Output_plots/Figure2c_signatures_fibros_p40.tiff", width=10, height=6, res=300, units='cm')

(ggplot(aes(y=value, x=factor(Sample,levels = samples_levels), fill = factor(variable, levels = sigs_levels)), data = result_df)
  + geom_bar(stat = "identity")
  + scale_y_continuous(labels = scales::percent) 
  + facet_grid(cols = vars(result_df$Type), scale = 'free_x',space = "free_x")
  + theme_bw() 
  + scale_fill_manual(values = c("darkred","bisque1","bisque3","darkgrey","beige","darkblue"),labels = c("SBS10c\n(POLD1 deficiency)", "SBS36", "SBS37", "SBS45\n(possible artefacts)", "SBS93", "SBS5\n(clock-like)"))
  + xlab("")
  + ylab("")
  + theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=7))
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  + theme(strip.text = element_text(size = 10))
  + theme(legend.key.height = unit(0.5, "cm"),legend.margin=margin(t = 0, unit='cm'),legend.title=element_blank())
)

dev.off()


#Figure 2d

df = read.table("Input_files/Figure2BC_fibros_mutrate_data.txt", header=T)

df_L474P = df[df$Status == "L474P",]
df_wt = df[df$Status == "wt",]

N_mut_wt = mean(df_wt$N_mut_all)
N_SBS5_wt = mean(df_wt$SBS5)
N_SBS10c_wt = mean(df_wt$SBS10c)

df_L474P_SBS5 = df_L474P[,c(1,3)]
colnames(df_L474P_SBS5) = c("Sample", "N_mut")
df_L474P_SBS5$Type = "SBS5"
df_L474P_SBS5$wt_value = N_SBS5_wt

df_L474P_SBS10c = df_L474P[,c(1,4)]
colnames(df_L474P_SBS10c) = c("Sample", "N_mut")
df_L474P_SBS10c$Type = "SBS10c"
df_L474P_SBS10c$wt_value = N_SBS10c_wt

result = rbind(df_L474P_SBS5,df_L474P_SBS10c)

samples_order = c('IV.4', 'III.4', 'III.2', 'IV.2')
sig_levels = c("SBS5", "SBS10c")

tiff(filename="Output_plots/Figure2d_fibros_mutrate_by_signature.tiff", width=6.5, height=6, res=300, units='cm')
(ggplot(result, aes(x=factor(Sample,levels=samples_order), y=N_mut, fill=Type))
  + geom_bar(stat="identity", alpha=0.8)
  + theme_bw()
  + facet_grid(cols = vars(factor(result$Type, levels = sig_levels)), scales = "free_y")
  + theme( panel.border = element_blank())
  + theme(panel.spacing.x = unit(1.5,"line"))
  + scale_fill_manual(values = c('darkblue','darkred'))
  + xlab("")
  + ylab("")
  + geom_hline(yintercept = result$wt_value, linetype = "dashed", col="black", size=1.2)
  + theme(legend.position = 'none',strip.text = element_text(size = 10))
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8)))
dev.off()

#Figures 2e-2g

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

#2f
tiff(filename="Output_plots/Figure2f_fibros_p40_pca.tiff", width=8, height=6, res=300, units='cm')
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

#2g

weights_pc_by_type<-as.data.frame(res.pca_by_type$rotation)
weights_pc_by_type$Mut_type = rownames(weights_pc_by_type)
weights_pc_by_type$Mutation = substr(weights_pc_by_type$Mut_type,3,5)
weights_pc_by_type$minus_PC1 = -weights_pc_by_type$PC1

muttype_level = c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T" )
x_labels = c("A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T","A.A","A.C","A.G","A.T","C.A","C.C","C.G","C.T","G.A","G.C","G.G","G.T","T.A","T.C","T.G","T.T" )


tiff(filename="Output_plots/Figure2e_PC1_spectrum.tiff", width=18, height=6, res=300, units='cm')
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

#2g

Tissue<-c("Fibros","Fibros","Fibros","Fibros_wt","Fibros","Fibros_wt")
result_pc_by_type$Tissue = Tissue

#download data for crypts and project to pc space
vcf_files_colon_PolD <- list.files(path = './Input_files/crypts/', pattern="*.vcf", full.names = T)
vcfs_colon_PolD <- read_vcfs_as_granges(vcf_files_colon_PolD, vcf_files_colon_PolD, ref_genome,  type = "snv")
mut_mat_colon_PolD <- mut_matrix(vcf_list = vcfs_colon_PolD, ref_genome = ref_genome)
mut_mat_freq_by_type_colon_PolD = create_mut_mat_freq_by_type(mut_mat_colon_PolD)
colon_PolD_pca_by_type<-as.data.frame(predict(res.pca_by_type, mut_mat_freq_by_type_colon_PolD))
colon_PolD_pca_by_type$Tissue = "colon"
colon_PolD_pca_by_type$Pol = "PolD"
colon_PolD_pca_by_type$PolD_status = c(rep("S478N", 35), rep("L474P", 8), rep("D316N",9))

#download data for sperm and project to pc space
vcf_files_duplex_PolD <- list.files(path = './Input_files/blood_sperm/', pattern="*.vcf", full.names = T)
vcfs_duplex_PolD <- read_vcfs_as_granges(vcf_files_duplex_PolD, vcf_files_duplex_PolD, ref_genome,  type = "snv")
mut_mat_duplex_PolD <- mut_matrix(vcf_list = vcfs_duplex_PolD, ref_genome = ref_genome)
mut_mat_freq_by_type_duplex_PolD = create_mut_mat_freq_by_type(mut_mat_duplex_PolD)
duplex_PolD_pca_by_type<-as.data.frame(predict(res.pca_by_type, mut_mat_freq_by_type_duplex_PolD))
duplex_PolD_pca_by_type$Tissue = c("sperm","blood")
duplex_PolD_pca_by_type$Pol = "PolD"
duplex_PolD_pca_by_type$PolD_status = "S478N"


#merge results for fibros, crypts and sperm in one dataframe
result_by_type = rbind(result_pc_by_type, colon_PolD_pca_by_type)
result_by_type = rbind(result_by_type, duplex_PolD_pca_by_type)

x_levels = c("Fibro_wt", "Fibro_PolD", "colon","blood", "sperm")
tissue_levels = c("Fibros_wt", "Fibros", "colon", "blood","sperm")
mutation_levels = c("wt", "L474P","D316N", "S478N")

tiff(filename="Output_plots/Figure2g_other_tissues_PC1.tiff", width=8, height=6, res=300, units='cm')
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
  + scale_x_discrete(labels=c("Fibros_wt" = "Fibros\nwt", "Fibros" = "Fibros\nL474P","colon" = "colon", "blood"="blood","sperm"="sperm"))
  + scale_color_discrete(name = "Mutation"))
dev.off()

