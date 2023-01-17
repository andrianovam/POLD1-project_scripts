library(ggridges)
library(viridis)
library(ggtext)
library(this.path)
library(MutationalPatterns)
library(gridExtra)
library(ggplot2)
library(lsa)
library(sigfit)
library(reshape2)
library(ggh4x)
library(reshape2)
library(RColorBrewer)
library(Rcpp)
library(ggpmisc)
library(grid)

setwd(dirname(this.path()))

#Figure 4a
df = read.table("Input_files/Figure4A_hyper_ultra_mutability.txt", header=T, sep="\t")
df = df[df$Status!="PolD_MSI-H-nonfunctional",]
df = df[df$Status!="PolD_PolE",]
status_levels = c("MSS", "PolE_MSS", "MSI", "PolD_MSI-H", "PolE_MSI-H", "POLE", "POLD1")
mut_per_mb_413711TT	= 226.0333333
mut_per_mb_second_tumour	= 96.4
mut_per_mb_hypermutable_polyps_stratton = mean(c(241677/3000,236160/3000, 233453/3000))


tiff(filename="Output_plots/Figure4a_Mutability_of_different_ucec_tcga_wes.tiff", width=12, height=7, res=300, units='cm')

(ggplot(df, aes(y=factor(Status,levels = status_levels), x=Mut_per_mb,fill=Tissue))
  + geom_density_ridges(alpha=0.7)+scale_x_continuous(trans='log',breaks=c(10,100))
  + theme_bw() 
  + scale_fill_brewer(palette = "Dark2")
  + geom_vline(xintercept = 10,linetype='dashed', col="darkgrey")
  + geom_vline(xintercept = 100,linetype='dashed', col="darkgrey")
  + annotate("point", x = mut_per_mb_413711TT, y = 7, col="#7570B3", size=2)
  + annotate('point', x = mut_per_mb_second_tumour, y = 7, col="#7570B3", size=2)
  + annotate('point', x = mut_per_mb_hypermutable_polyps_stratton, y = 7, col='#D95F02', size=2)
  + ylab("")
  + xlab("Mutations per megabase")
  + annotate("text", x = 1, y = 8.5,label = "Nonmutable", size=2.9)
  + annotate("text", x = 33, y = 8.5, label = "Hypermutable", size=2.9)
  + annotate("text", x = 500, y = 8.5, label = "Ultrahypermutable", size=2.9)
  + annotate("text", x = 100, y = 7.8, label = "S478N", size=2.5, col='#D95F02', fontface=2)
  + annotate("text", x = 140, y = 7.35, label = "D316H", size=2.5, col='#7570B3', fontface=2)
  + annotate("text", x = 370, y = 7.35, label = "L474P", size=2.5, col='#7570B3', fontface=2)
  + theme(legend.position = 'top')
  + theme(legend.key.size = unit(0.4, "cm"))
  + theme(legend.margin =margin(r=10,l=5,t=0,b=0))
  + theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=8),legend.title=element_text(size=10))
  + theme(axis.text.y = element_markdown())
  + scale_y_discrete(expand = expansion(mult = c(0, 0.3)), labels=c("POLD1" = "*POLD1exo-*", "POLE" = "*POLEexo-*", "PolE_MSI-H" = "MMRd *POLEexo-*", "PolD_MSI-H" = "MMRd *POLD1exo-*", "MSI" = "MMRd", "PolE_MSS" = "*POLEexo-*", "MSS"="MMR proficiency<br>proofread proficiency")))

dev.off()


#Figure 4b

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

vcf_files <- c("Input_files/Tumors/IV_1_T.target_more_0.15.vcf","Input_files/Tumors/D316H_T_somatic_filtered_0.15_PASS.vcf")
vcf_files
sample_names <- c("POLD1_L4747P","POLD1_D316H")
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome,  type = "all")
summary(vcfs)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
type_occurrences

mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

tiff(filename="Output_plots/Figure4b_tumor_spectra.tiff", width=12, height=7, res=300, units='cm')
plot_96_profile(mut_mat, condensed = TRUE)
dev.off()


#Figure 4c
data("cosmic_signatures_v3.2")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

vcf_files <- list.files(path = "Input_files/crypts", pattern="*.vcf", full.names = T)
sample_names <- list.files(path = "Input_files/crypts", pattern="*.vcf")
vcf_files_tumors <- list.files(path = "Input_files/Tumors/", full.names = T)
sample_names_tumors <- c("IV.1_T","D316H_T")

vcf_files_all = c(vcf_files, vcf_files_tumors)
sample_names_all = c(sample_names, sample_names_tumors)


vcfs <- read_vcfs_as_granges(vcf_files_all, sample_names_all, ref_genome,  type = "snv")
summary(vcfs)

tissue <- c(rep("S478N", 21), rep("S478N_hypermutable",3), rep("S478N", 11), rep("L474P", 8), rep("D316N", 9), "L474P_hypermutable", "D316H_hypermutable")
mut_mat_crypts <- t(mut_matrix(vcf_list = vcfs, ref_genome = ref_genome))
mut_mat_crypts=as.data.frame(mut_mat_crypts)
mut_mat_crypts$Mutation = tissue
test = aggregate(.~Mutation,data = mut_mat_crypts,  sum)
rownames(test)=test$Mutation
test = test[,c(2:ncol(test))]
mcmc_samples_fit <- fit_signatures(counts = test, 
                                   signatures = cosmic_signatures_v3.2[c(1,5,15,16),],
                                   iter = 5000, 
                                   warmup = 1000, 
                                   chains = 1, 
                                   seed = 1756)

exposures <- retrieve_pars(mcmc_samples_fit, 
                           par = "exposures", 
                           hpd_prob = 0.90)
names(exposures)

exposures_mean = exposures$mean
exposures_mean$Mutation = rownames(exposures_mean)
exposures_melt = melt(exposures_mean)


mutation_levels = c("D316N", "L474P", "S478N", "S478N_hypermutable", "L474P_hypermutable","D316H_hypermutable")

tiff(filename="Output_plots/Figure4?_crypts_tumors_somatic_signatures.tiff", width=10, height=6, res=300, units='cm')
(ggplot(exposures_melt, aes(x=factor(Mutation, mutation_levels), y=value, fill=variable))
  +geom_bar(stat="identity",width=0.7)
  +theme_bw()
  +ylab("Mutation fraction")
  +scale_fill_manual(values = c("lightcoral","indianred4","plum2","purple4"),name="")
  +scale_x_discrete(labels=c("D316N" = "D316N\nintestinal\ncrypts", "L474P" = "L474P\nintestinal\ncrypts","S478N" = "S478N\nintestinal\ncrypts", "S478N_hypermutable"="S478N\nAdenomas", "L474P_hypermutable"="L474P\nCancer","D316H_hypermutable"="D316H\nCancer"))
  +xlab("")
  +theme(axis.title.x=element_blank())
  +theme(legend.position = 'top')
  +theme(legend.key.size = unit(0.4, "cm"))
  +theme(legend.margin =margin(r=10,l=5,t=5,b=0))
  +theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
  +theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=8),legend.title=element_text(size=10)))
dev.off()


#Figure 4d

df = read.table("Input_files/413701TT_chr19_germline_vars.txt")
df = df[,c(1:5)]
df$vaf1 = df$V4/df$V3
df$vaf2 = df$V4/df$V3
colnames(df)<-c("Chr","Pos", "Cov","Ref","Alt","Vaf1","Vaf2")
df$Cov_log = log(df$Cov)
df_vaf1 = df[,c(1,2,6)]
colnames(df_vaf1)<-c("Chr","Pos","Value")
df_vaf2 = df[,c(1,2,7)]
colnames(df_vaf2)<-c("Chr","Pos","Value")
df_for_plot_vafs = rbind(df_vaf1, df_vaf2)
df_for_plot_vafs$Type="vaf"
df_cov = df[,c(1,2,8)]
colnames(df_cov)<-c("Chr","Pos", "Value")
df_cov$Type="Cov,log"
df_for_plot = rbind(df_for_plot_vafs, df_cov)

mut_vaf = df[df$Pos==50909701,]$Vaf
mut_pos = 50909701


tiff(filename="Output_plots/Figure4d_413701TT_chr19_allelic_disbalance.tiff", width=9, height=6, res=300, units='cm')
(ggplot(df_for_plot, aes(x=Pos, y=Value))
  +geom_point(size=0.5)
  +geom_vline(xintercept = mut_pos,col="red", size=1)
  +theme_bw()
  +facet_grid(rows = vars(factor(Type, levels=c("vaf","Cov,log"))), scale="free_y",switch = "y")
  +force_panelsizes(rows = c(1, 0.6))+ylab("") 
  +theme(panel.spacing.y = unit(0, "mm"),strip.background = element_rect(size = 0.5))
  +ggtitle("*POLD1* L474P Ultrahypermutable colon cancer")
  +theme(legend.position='none')
  +theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=8),legend.title=element_text(size=10),plot.title = element_markdown(size=8))
  +theme(plot.margin = margin(0.2, 0.5, 0.5, 0.5, "cm"))
  +xlab("Chr19 position"))
dev.off()


#Figure 4e

df = read.table("Input_files/100kb_windows_pold_S478N_context_CtoA_TCT.txt")
colnames(df) = c("chr", "start", "end", "id", "N_mut")
df$Type = "wt/S478N"
df1 = read.table("Input_files/100kb_windows_pold_S478N_hypermutable_context_CtoA_TCT.txt")
colnames(df1) = c("chr", "start", "end", "id", "N_mut")
df1$Type = "S478N/S478N"
df_sites = read.table("Input_files/100kb_windows_TCT_AGA_sites.txt")
colnames(df_sites) = c("id", "sites")
df = merge(df, df_sites, by="id")
df$sites_all = df$sites * 32
df1 = merge(df1, df_sites, by="id")
df1$sites_all = df$sites * 1
result = rbind(df, df1)

df_rt = read.table("Input_files/HeLa-S3_hg19_100kb_windows.bed")
colnames(df_rt)=c("id","V2","V3","V4","V5","RT")
quantiles_rt = quantile(df_rt$RT, probs = seq(0,1,0.02))
quantile_bin <-function(x){
  bin=0
  for (i in seq(1,50)){
    if (x > quantiles_rt[i] & x <= quantiles_rt[i+1]){
      bin = i
    }
  }
  return(bin)
} 

df_rt$rt_bin<-as.numeric(lapply(df_rt$RT, quantile_bin))

result = merge(result, df_rt[,c(1,7)], by="id")

x = aggregate(result[,c(5,8)], list(result$Type,result$rt_bin), FUN=sum)
colnames(x) = c("Mut_type","rt_bin","N_mut","Sites")
x$Nonmut = x$Sites-x$N_mut
model <- glm(data = x,cbind(N_mut, Nonmut) ~ Mut_type+rt_bin+Mut_type*rt_bin, family="binomial")

homo = x[x$Mut_type=='S478N/S478N',]
homo = homo[homo$rt_bin >3,]
homo$mutrate = homo$N_mut/homo$Sites
homo$mutrate_norm = homo$mutrate/min(homo$mutrate)

hetero = x[x$Mut_type=='wt/S478N',]
hetero = hetero[hetero$rt_bin >3,]
hetero$mutrate = hetero$N_mut/hetero$Sites
hetero$mutrate_norm = hetero$mutrate/min(hetero$mutrate)

y = rbind(homo,hetero)

tiff(filename="Output_plots/Figure4e_mutrate_vs_rt_homo_hetero.tiff", width=8, height=7, res=300, units='cm')
(ggplot(y, aes(x=rt_bin, y=mutrate_norm, group=Mut_type,col=Mut_type))
  + geom_point()
  + geom_smooth(method="lm") 
  + theme_bw()
  + scale_color_manual(values=c("purple4","plum3"), name="POLD1")
  + theme(legend.position = 'bottom')
  + xlab("Late <--- Replication timing ---> Early") 
  + ylab("TCT>TAT mutation rate")
  + theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=8),legend.title=element_text(size=10))
  + stat_poly_eq(aes(label = after_stat(eq.label)),label.x = 0.95, vstep=0.15, size=4)
  +  annotate(geom="text", x=40, y=4, label=paste("p-val","<2e-16", sep=""), size=3.5)
  )
dev.off()


#Figure 4f

df = read.table("Input_files/Figure4F_tcga_ucec_mutrate.txt", header=T, sep="\t")

df_aggregate = do.call(data.frame, aggregate(df$nSNV_from_vcf, by=list(df$Status_MMR, df$Status_polymerase), FUN = function(x) c(mean = mean(x), median = median(x), sd = sd(x))))
df_aggregate=rbind(df_aggregate, c("MSS", "POLD1", 0,0,0))
colnames(df_aggregate) = c("Status_MMR", "Status_polymerase", "Mean_mutrate","Median_mutrate", "SD_mutrate")
df_aggregate$Median_mutrate=as.numeric(df_aggregate$Median_mutrate)
df_aggregate$Mean_mutrate=as.numeric(df_aggregate$Mean_mutrate)
df_aggregate$SD_mutrate=as.numeric(df_aggregate$SD_mutrate)
df_aggregate$Median_mutrate = as.numeric(df_aggregate$Median_mutrate)
head(df_aggregate)


tiff(filename="Output_plots/Figure4F_tcga_ucec_mutrate.tiff", width=8, height=7, res=300, units='cm')
(ggplot(df_aggregate, aes(x=Status_polymerase, y=(Median_mutrate/30),fill=Status_MMR))
  +geom_bar(stat="identity", position = 'dodge')
  +theme_bw()
  +theme(legend.position="bottom")
  +ylab("Median mutation rate per Mb")
  +xlab("")
  +theme(axis.text.x = element_markdown())
  +scale_x_discrete(labels=c("POLD1" = "*POLD1*<br>mutated", "POLE" = "*POLE*<br>mutated","wt" = "wt<br>polymerases"))
  +scale_fill_manual(values = c("coral3","cadetblue3"),name="")
  +theme(axis.title.x=element_blank())
  +theme(legend.key.size = unit(0.4, "cm"))
  +theme(legend.margin =margin(r=10,l=5,t=0,b=0))
  +theme(plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"))
  +theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=8),legend.title=element_text(size=10)))
dev.off()

