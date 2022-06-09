library("MutationalPatterns")
library("gridExtra")
library("ggplot2")
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)

setwd(dirname(this.path()))

percent_of_pold_samples = 0 
n_simulations_per_family = 1

#### create spectrum of de novo mutations in offspring of fathers with POLD1 L474P variant###
setwd("../Figures/Input_files/de_novo")
vcf_files_germline <- c("IV.4_SNP_mv20.vcf","IV.5_SNP_mv20.vcf","IV.6_SNP_mv20.vcf","IV.8_SNP_mv20.vcf")
vcf_files_germline
sample_names_germline <- c("IV.4", "IV.5","IV.6","IV.8")
vcfs_germline <- read_vcfs_as_granges(vcf_files_germline, sample_names_germline, ref_genome,  type = "all")
summary(vcfs_germline)
pold_spectrum <- mut_matrix(vcf_list = vcfs_germline, ref_genome = ref_genome)
pold_spectrum = as.data.frame(pold_spectrum)
pold_spectrum$sum = rowSums(pold_spectrum)
pold_spectrum$Freq = pold_spectrum$sum/(sum(pold_spectrum$sum))
pold_spectrum$Mutation = rownames(pold_spectrum)


### create spectrum of de novo mutations in public dataset of trios ###

setwd("..")
df = read.table("public_trio_data_matrix.txt")
df$N_mut = rowSums(df)
wt_spectrum = data.frame(Mutation = rownames(df), N_mut = df$N_mut)
wt_spectrum$Freq = wt_spectrum$N_mut/sum(wt_spectrum$N_mut)
head(wt_spectrum)

families_id = colnames(df)
n_mutations_per_family = colSums(df)
  
### prepare simulations ###
result<-NULL
all_simulations = length(n_mutations_per_family) * n_simulations_per_family
n_simulations_done = 0
n_simulations_norm = round(all_simulations * (1-percent_of_pold_samples), 0)

for (i in c(1:length(n_mutations_per_family))){
  print(i)
  number_of_mutations = n_mutations_per_family[i]
  for (j in c(1:n_simulations_per_family)){    #number of simulations per family
    n_simulations_done = n_simulations_done + 1
    print(j)
    quasi_id = paste(families_id[i], as.character(j), sep="_")
    if (n_simulations_done <= n_simulations_norm){
    subsample = sample(wt_spectrum$Mutation, size=number_of_mutations, replace=TRUE, prob=wt_spectrum$Freq)} else{
      subsample = sample(pold_spectrum$Mutation, size=number_of_mutations, replace=TRUE, prob=pold_spectrum$Freq)    
    }
    df_cur = data.frame(mut = subsample, family_id = quasi_id)
    result=rbind(result, df_cur)
  }
}

result$Mut_type = substr(result$mut,3,5)
head(result)

###convert result table to table for pca analysis ###

families = unique(result$family_id)

mut_mat_freq_by_type_simulations<-NULL
ind = 0
for (family in families){
  print(family)
  ind = ind+1
  print(ind)
  df_family = result[result$family_id == family,]
  family_result = NULL
  possible_mutations = c("A[C>A]A","A[C>A]C","A[C>A]G","A[C>A]T","C[C>A]A","C[C>A]C","C[C>A]G","C[C>A]T","G[C>A]A","G[C>A]C","G[C>A]G","G[C>A]T","T[C>A]A","T[C>A]C","T[C>A]G","T[C>A]T","A[C>G]A","A[C>G]C","A[C>G]G","A[C>G]T","C[C>G]A","C[C>G]C","C[C>G]G","C[C>G]T","G[C>G]A","G[C>G]C","G[C>G]G","G[C>G]T","T[C>G]A","T[C>G]C","T[C>G]G","T[C>G]T","A[C>T]A","A[C>T]C","A[C>T]G","A[C>T]T","C[C>T]A","C[C>T]C","C[C>T]G","C[C>T]T","G[C>T]A","G[C>T]C","G[C>T]G","G[C>T]T","T[C>T]A","T[C>T]C","T[C>T]G","T[C>T]T","A[T>A]A","A[T>A]C","A[T>A]G","A[T>A]T","C[T>A]A","C[T>A]C","C[T>A]G","C[T>A]T","G[T>A]A","G[T>A]C","G[T>A]G","G[T>A]T","T[T>A]A","T[T>A]C","T[T>A]G","T[T>A]T","A[T>C]A","A[T>C]C","A[T>C]G","A[T>C]T","C[T>C]A","C[T>C]C","C[T>C]G","C[T>C]T","G[T>C]A","G[T>C]C","G[T>C]G","G[T>C]T","T[T>C]A","T[T>C]C","T[T>C]G","T[T>C]T","A[T>G]A","A[T>G]C","A[T>G]G","A[T>G]T","C[T>G]A","C[T>G]C","C[T>G]G","C[T>G]T","G[T>G]A","G[T>G]C","G[T>G]G","G[T>G]T","T[T>G]A","T[T>G]C","T[T>G]G","T[T>G]T" )
  for (i in possible_mutations){
    cur_mut = substr(i,3,5)
    sum_cur_mut = sum(df_family$Mut_type == cur_mut) 
    if (sum_cur_mut !=0){
      context_freq = sum(df_family$mut == i)/sum_cur_mut} else{
        context_freq = 0
      }
    context_result = c(i, context_freq)
    family_result = rbind(family_result, context_result)
  }
  family_result=as.data.frame(family_result)
  family_result[,2] = as.numeric(family_result[,2])
  colnames(family_result)<-c("mut", family)
  if (is.null(mut_mat_freq_by_type_simulations)==TRUE){
    mut_mat_freq_by_type_simulations<-family_result}else
    {mut_mat_freq_by_type_simulations<-merge(mut_mat_freq_by_type_simulations,family_result, by="mut")}
}


rownames(mut_mat_freq_by_type_simulations)<-mut_mat_freq_by_type_simulations$mut
mut_mat_freq_by_type_simulations<-mut_mat_freq_by_type_simulations[2:ncol(mut_mat_freq_by_type_simulations)]
mut_mat_freq_by_type_simulations[1:5,1:5]


#write.table(mut_mat_freq_by_type_simulations, file = paste("matrix_for_pca_simulations_prop_pold_", percent_of_pold_samples,".txt",sep=""), quote = F)
