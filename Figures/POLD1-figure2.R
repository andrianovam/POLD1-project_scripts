library(ggplot2)
library(ggtext)
library(dplyr)
library("MutationalPatterns")
library("this.path")
library(cowplot)
library(gridExtra)

setwd(dirname(this.path()))

#Figure 2a

df = read.table("Input_files/Figure2A_signatures_fibros_p0.txt", header=T)
result_df = melt(df)

sigs_levels = c("SBS58","SBS10c","SBS1", "SBS5", "SBS7a", "SBS7b", "SBS7d")
samples_levels = c('III.2', 'III.4','IV.1','IV.2', 'IV.4','IV.6', 'IV.3', 'IV.5')

tiff(filename="Output_plots/Figure2a_signatures_fibros_p0.tiff", width=9, height=6, res=300, units='cm')

(ggplot(aes(y=value, x=factor(Samples,levels = samples_levels), fill = factor(variable, levels = sigs_levels)), data = result_df)
  + geom_bar(stat = "identity")
  + facet_grid(cols = vars(result_df$Status), scale = 'free_x',space = "free_x")
  + theme_bw() 
  + scale_fill_manual(values = c("grey", "darkred", "blue", "darkblue","darkseagreen1", "darkseagreen2", "darkseagreen3"), name="")
  + xlab("")
  + ylab("")
  + theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=8))
  + theme(legend.key.size = unit(0.4, "cm"))
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  + theme(strip.text = element_text(size = 10)))

dev.off()

#Figure 2b

df = read.table("Input_files/Figure2BC_fibros_mutrate_data.txt", header=T)

samples_levels = c("IV.5","IV.3", "IV.2", "III.2", "III.4", "IV.4")

tiff(filename="Output_plots/Figure2b_mutrate_fibros.tiff", width=7.5, height=5, res=300, units='cm')
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

tiff(filename="Output_plots/Figure2c_signatures_fibros_p40.tiff", width=9, height=6, res=300, units='cm')

(ggplot(aes(y=value, x=factor(Sample,levels = samples_levels), fill = factor(variable, levels = sigs_levels)), data = result_df)
  + geom_bar(stat = "identity")
  + scale_y_continuous(labels = scales::percent) 
  + facet_grid(cols = vars(result_df$Type), scale = 'free_x',space = "free_x")
  + theme_bw() 
  + scale_fill_manual(values = c("darkred","bisque1","bisque3","grey","beige","darkblue"), name="")
  + xlab("")
  + ylab("")
  + theme(axis.text=element_text(size=8),axis.title=element_text(size=10),legend.text=element_text(size=8))
  + theme(legend.key.size = unit(0.4, "cm"))
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  + theme(strip.text = element_text(size = 10)))

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

tiff(filename="Output_plots/Figure2d_fibros_mutrate_by_signature.tiff", width=9, height=5, res=300, units='cm')
(ggplot(result, aes(x=factor(Sample,levels=samples_order), y=N_mut, fill=Type))
  + geom_bar(stat="identity", alpha=0.7)
  + geom_hline(yintercept = result$wt_value, linetype = "dashed", col="black", size=1.5)
  + theme_bw()
  + facet_grid(cols = vars(factor(result$Type, levels = sig_levels)), scales = "free_y")
  + theme( panel.border = element_blank())
  + theme(panel.spacing.x = unit(1.5,"line"))
  + scale_fill_manual(values = c('darkblue','darkred'))
  + xlab("")
  + ylab("")
  + theme(legend.position = 'none',strip.text = element_text(size = 10))
  + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  + theme(axis.text=element_text(size=9),axis.title=element_text(size=10),legend.text=element_text(size=8)))
dev.off()

