library(ggplot2)
library(ggtext)
library(dplyr)
library("MutationalPatterns")
library("this.path")
library(cowplot)
library(gridExtra)

setwd(dirname(this.path()))

#Figure 3a

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

tiff(filename="Output_plots/Figure3a_germline_counts.tiff", width=8, height=8, res=300, units='cm')

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
  + scale_x_discrete(labels=c("Father" = "Father\nL474P", "Mother" = "Mother\nL474P", "Both wt" = "Both wt\n(our data)", "Published trios" = "Both wt\n(published)"))
  + stat_summary(data = data, aes(x=factor(POLD,levels = levels_pold),y=mut), fun.data=mean_sdl, geom="pointrange",size=0.1, col='black')
  + theme(legend.position='none')
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8))
)

dev.off()

#Figure 3b
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

tiff(filename="Output_plots/Figure3C_germline_4contexts_prop.tiff", width=8, height=8, res=300, units='cm')
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
  + scale_x_discrete(labels=c("Father" = "Father\nL474P", "Mother" = "Mother\nL474P", "Both wt" = "Both wt\n(our data)", "Published trios" = "Both wt\n(published)"))
  + stat_summary(data = data1, aes(x=factor(POLD,levels = levels_pold),y=four_contexts_proportion), fun.data=mean_sdl, geom="pointrange", color="black",size=0.2)
  + theme(legend.position='none')
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8))
)

dev.off()

#Figure 3c
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

tiff(filename="Output_plots/Figure3c_germline_PC1_by_type.tiff", width=8, height=8, res=300, units='cm')


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
  + scale_x_discrete(labels=c("Father" = "Father\nL474P", "Mother" = "Mother\nL474P", "Both wt" = "Both wt\n(our data)", "Published trios" = "Both wt\n(published)"))
  + stat_summary(data = data1, aes(x=factor(Status,levels = levels_pold),y=-PC1), fun.data=mean_sdl, geom="pointrange", color="black",size=0.2)
  + theme(legend.position='none')
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8))
)

dev.off()

#Figure 3d-f
df = read.table("Input_files/public_trio_data.txt")
head(df)
colnames(df) = c("chr","coord","mut","family_id","age","dataset")
count_mut_per_family = aggregate(mut ~ family_id, data = df[,c(3,4)], FUN = length)
head(count_mut_per_family)
our_data = data.frame(sample = c("IV.1","IV.2","IV.3","IV.4","IV.5","IV.6","IV.8","IV.9"), mut=c(64,62,71,61,64,88,24,50), PolD_carrier = c("mother POLD1 L474P","mother POLD1 L474P","mother POLD1 L474P","father POLD1 L474P","father POLD1 L474P","father POLD1 L474P","father POLD1 L474P", "wt parents"), nudge_x = c(3,0,2,-4,3,0,0,0),nudge_y = c(0,0,0,0,30,0,0,0))

ks_pval_fathers = round(ks.test(count_mut_per_family$mut, our_data[our_data$PolD_carrier=="father POLD1 L474P",]$mut)$p.value,4)
ks_pval_mothers = round(ks.test(count_mut_per_family$mut, our_data[our_data$PolD_carrier=="mother POLD1 L474P",]$mut)$p.value,4)

#plot number of mutations
p1<-(ggplot(count_mut_per_family, aes(x=mut)) 
     + geom_histogram(fill="grey",bins=5000,binwidth = 1,alpha=0.5) 
     + theme_bw() 
     + geom_point(data=our_data, aes(x = mut, y=0, col=PolD_carrier),size=2)
     + xlab("N mutations per genome")
     + ylab("")
     + scale_color_manual(values=c("indianred3","khaki3","black"))
     + geom_text(data = our_data,aes(label=sample, x=mut,y=30), size=3,show.legend=FALSE, nudge_x=our_data$nudge_x, nudge_y = our_data$nudge_y)
     + xlim(0,130)
     + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none')
     + annotate("text", x=120, y=170, label= paste("p-val", ks_pval_fathers, sep="="), color = 'indianred3')
     + annotate("text", x=120, y=130, label= paste("p-val", ks_pval_mothers, sep="="), color = 'khaki3')
     + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8)))

#plot PC1
normal_germline_pca = read.table("Input_files/public_germline_pca_by_type.txt", header=T)
our_germline_pca = read.table("Input_files/our_germline_pca_by_type.txt", header=T)
our_germline_pca$PolD_carrier = c("Father","Father","Father","Mother","Father","Mother","Mother", "wt parents")
our_germline_pca$col = c("red","red","red","khaki3","red","khaki3","khaki3", "black")
our_germline_pca$Sample = rownames(our_germline_pca)
our_germline_pca$nudge = c(-0.03,0,0,0,0.03,-0.03,0.03,0)
our_germline_pca$Sample_name = c('IV.8','IV.6','IV.5','IV.3','IV.4','IV.2','IV.1','IV.9')

ks_pval_PC1_fathers = round(ks.test(normal_germline_pca$PC1, our_germline_pca[our_germline_pca$PolD_carrier=="Father",]$PC1)$p.value,4)
ks_pval_PC1_mothers = round(ks.test(normal_germline_pca$PC1, our_germline_pca[our_germline_pca$PolD_carrier=="Mother",]$PC1)$p.value,4)


p2<-(ggplot(normal_germline_pca, aes(x=-PC1)) 
     + geom_histogram(fill="grey",bins=100,alpha=0.5) 
     + theme_bw() 
     + geom_point(data=our_germline_pca, aes(x = -PC1, y = 0,col=PolD_carrier),size=2)
     + geom_vline(xintercept = quantile(-normal_germline_pca$PC1, probs = seq(0,1,0.01))[91], col="black", size=0.8)
     + scale_color_manual(values=c("indianred3","khaki3","black"))
     + geom_text(data = our_germline_pca,aes(label=Sample_name, x=-PC2,y=30), size=3,show.legend=FALSE, nudge_x=our_germline_pca$nudge)
     + geom_text(label="90th percentile", x=0.13,y=200, size=3,show.legend=FALSE)
     + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'none')
     + ylab("")
     + annotate("text", x=0.45, y=190, label= paste("p-val=", ks_pval_PC1_fathers, sep=""),color = 'indianred3')
     + annotate("text", x=0.45, y=150, label= paste("p-val=", ks_pval_PC1_mothers, sep=""),color='khaki3')
     + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8)))


#plot cosine with SBS10c
our_germline_SBS10c_cosine = read.table("Input_files/our_germline_SBS10c_cosine.txt", header=T)
normal_germline_SBS10c_cosine = read.table("Input_files/public_germline_SBS10c_cosine.txt", header=T)
our_germline_SBS10c_cosine$Carrier = c("Father L474P","Father L474P","Father L474P","Mother L474P","Father L474P","Mother L474P","Mother L474P","wt parents")
our_germline_SBS10c_cosine$Sample = rownames(our_germline_SBS10c_cosine)
our_germline_SBS10c_cosine$nudge = c(-0.01,-0.02,0,0.01,0,0.01,-0.006,0)
our_germline_SBS10c_cosine$Sample_name = c('IV.8','IV.6','IV.5','IV.3','IV.4','IV.2','IV.1','IV.9')

ks_pval_SBS10c_cosine_fathers = round(ks.test(normal_germline_SBS10c_cosine$coef_SBS10c, our_germline_SBS10c_cosine[our_germline_SBS10c_cosine$Carrier=="Father L474P",]$coef_SBS10c)$p.value,4)
ks_pval_SBS10c_cosine_mothers = round(ks.test(normal_germline_SBS10c_cosine$coef_SBS10c, our_germline_SBS10c_cosine[our_germline_SBS10c_cosine$Carrier=="Mother L474P",]$coef_SBS10c)$p.value,4)


p4<-(ggplot(normal_germline_SBS10c_cosine, aes(x=coef_SBS10c)) 
     + geom_histogram(fill="grey",bins=100,alpha=0.5) 
     + theme_bw() 
     + geom_point(data=our_germline_SBS10c_cosine, aes(x = coef_SBS10c, y = 0,col=Carrier),size=2)
     + geom_vline(xintercept = quantile(normal_germline_SBS10c_cosine$coef_SBS10c, probs = seq(0,1,0.01))[91], col="black", size=0.8)
     + scale_color_manual(values=c("indianred3","khaki3","black"),name = "Parents *POLD1* status")
     + geom_text(data = our_germline_SBS10c_cosine,aes(label=Sample_name, x=coef_SBS10c,y=50), size=3,show.legend=FALSE, nudge_x=our_germline_SBS10c_cosine$nudge)
     + geom_text(label="90th percentile", x=0.38,y=290, size=3,show.legend=FALSE)
     + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = 'bottom')
     + theme(legend.title = element_markdown())
     + ylab("")
     + xlab("cosine similarity with SBS10c")
     + annotate("text", x=0.6, y=250, label= paste("p-val", ks_pval_SBS10c_cosine_fathers, sep="="), color = 'indianred3')
     + annotate("text", x=0.6, y=200, label= paste("p-val", ks_pval_SBS10c_cosine_mothers, sep="="), color = 'khaki3')
     + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=10)))

tiff(filename="Output_plots/Figure3d-e_germline_nmut_PC1_SBS10_by_type.tiff", width=17, height=12, res=300, units='cm')
plot_grid(p1, p2, p4, align = "v", nrow = 3, rel_heights = c(0.3, 0.3, 0.4))
dev.off()

#Figure 3g

df_trios = read.table("Input_files/public_germline_pca_by_type.txt", header=T)
df_trios$type = "Data"
df_simulations = read.table("Input_files/simulations_pca_by_type_prop_pold_0.txt")
df_simulations$type = "Simulation wt trios"
result = rbind(df_trios, df_simulations)
nrow(result)
type_levels = c("Data", "Simulation wt trios")
p1<-(ggplot(result, aes(x=-PC1, fill=factor(type, levels = type_levels)))
     +geom_density(alpha = 0.5)
     +theme_bw()
     + ylim(0,3.8)
     +scale_fill_manual(values=c("blue", "darkgrey"), name=NULL)
     +theme(legend.position='top')
     +theme(legend.key.size = unit(0.4, "cm"))
     +theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8)))

qq.out <- qqplot(x=-df_simulations$PC1, y=-df_trios$PC1, plot.it=FALSE)
qq.out <- as.data.frame(qq.out)
xylim <- range( c(qq.out$x, qq.out$y) )
p2<-(ggplot(qq.out, aes( x= x, y = y)) 
     + geom_point(size=0.5)
     + geom_abline( intercept=0, slope=1)
     + coord_fixed(ratio = 1, xlim=xylim, ylim = xylim)
     + xlab("Simulation")
     + ylab("Data")
     + theme_bw()
     + theme(axis.text=element_text(size=7),axis.title=element_text(size=8),legend.text=element_text(size=6)))

jpeg(filename=paste("Output_plots/Figure3g_PC1_normal_trios_vs_simuations.tiff", sep=""), width=9, height=9, res=300, units='cm')
p1 + annotation_custom(ggplotGrob(p2), xmin = 0, xmax = 0.7, ymin = 1.5, ymax = 3.9)
dev.off()

#Figure 3h

df_trios = read.table("Input_files/simulations_pca_by_type_prop_pold_0.05.txt",header=T)
df_trios$type = "Simulation with POLD1 mutants"
df_simulations = read.table("Input_files/simulations_pca_by_type_prop_pold_0.txt")
df_simulations$type = "Simulation wt trios"
result = rbind(df_trios, df_simulations)
nrow(result)
type_levels = c("Simulation with POLD1 mutants", "Simulation wt trios")
p1<-(ggplot(result, aes(x=-PC1, fill=factor(type, levels = type_levels)))
     + geom_text(x=-0.35, y=3.8, label="5% POLD1 L474P samples", size=3)
     + geom_density(alpha = 0.5)
     + ylim(0,3.8)
     + theme_bw()
     + scale_fill_manual(values=c("blue", "darkgrey"), name=NULL)
     + theme(legend.position='top')
     + theme(legend.key.size = unit(0.4, "cm"))
     + geom_segment(x = 0.4, y=1, xend = 0.27, yend = 0.3, arrow = arrow(length = unit(0.1, "cm")),lineend = c('round'),linejoin = c('round'),size=1)
     + geom_text(label='distribution shift', x=0.4,y=1.1, size=3,show.legend=FALSE)
     + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8)))

qq.out <- qqplot(x=-df_simulations$PC1, y=-df_trios$PC1, plot.it=FALSE)
qq.out <- as.data.frame(qq.out)
xylim <- range( c(qq.out$x, qq.out$y) )
p2<-(ggplot(qq.out, aes( x= x, y = y)) 
     + geom_point(size=0.5)
     + geom_abline( intercept=0, slope=1)
     + coord_fixed(ratio = 1, xlim=xylim, ylim = xylim)
     + xlab("Simulation wt trios")
     + ylab("Simulation with POLD1 mutants")
     + theme(axis.title.y = element_markdown())
     + theme_bw()
     + theme(axis.text=element_text(size=7),axis.title=element_text(size=6),legend.text=element_text(size=6)))


#setwd("D:/Lab/cancer/Article_PolD_family/suppl_plots")
#jpeg(filename=paste("Suppl_Figure2_PC2_germline_normal_trios_vs_simuations.jpeg", sep=""), width=10, height=4, res=300, units='in')
#grid.arrange(p1,p2,nrow=1)
#dev.off()

setwd("D:/Lab/cancer/Article_PolD_family/plots")
jpeg(filename=paste("Figure3H_simuations_norm_vs_simulations_pold_pca_6000trios.tiff", sep=""), width=9, height=9, res=300, units='cm')
p1 + annotation_custom(ggplotGrob(p2), xmin = 0.05, xmax = 0.7, ymin = 1.15, ymax = 4)
dev.off()

