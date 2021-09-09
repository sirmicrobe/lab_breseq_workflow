#!/mmfs1/groups/HPC-Marshall/miniconda3/bin/Rscript --vanilla

#or if running locally would do #!/usr/bin/env Rscript --vanilla
#name: reformat gdtools output from breseq to combine all samples into one file
#built for output.gd files based on subtract and then annotate outputs from GDtools in breseq v0.35.7

#input position 1: full path to directory containing all breseq $i_subtractref_name.txt files
#input position 2: sample key in .csv format. must contain a column called "title" that matches the breseq number

library("tidyverse")


#command line inputs
args <- commandArgs(trailingOnly=TRUE)

# Save all filenames in directory ending in .txt as a variable temp
setwd(args[1])
temp = list.files(pattern="*.txt") 
info <- file.info(temp)
empty <- rownames(info[info$size <=1,])
if (length(empty)>=1) {
    temp <- temp[!temp == empty]
}
all_subinh <- lapply(temp,read.delim) #create a list with each object in list a data frame of mutation table

# headers we want to keep
headers <- c("aa_new_seq","aa_position","aa_ref_seq","codon_new_seq","codon_number","codon_position","codon_ref_seq","start","end","side_1_position","side_2_position","locus_tag","gene_name","gene_position","gene_product","frequency","mutation_category","new_seq","position","seq_id","size","snp_type","type","title")
#take a look at the first data frame in the list
df_all <- all_subinh[[1]]
df_all <- df_all[,headers] #select only headers we want to keep
head(df_all)
# convert all data frames in the list to one large data frame
all_df <- as.data.frame(data.table::rbindlist(all_subinh,fill=T),colClasses = c("character"))
all_df <- all_df[,headers] #keep headers we want
all_df <- all_df %>% filter(type!="UN") #remove UN at the end
all_df$mutation <- paste(all_df$aa_ref_seq,all_df$aa_position,all_df$aa_new_seq, sep=":") #create mutation column


#get rid of incompatible characters in file for R
system(paste("iconv -c -f utf-8 -t ascii",args[2], "> SampleKey_noutf.csv",sep=" "))

sample_key <- read.csv("SampleKey_noutf.csv",header=T,check.names=T)


all_df <- (merge(all_df,sample_key,by="title"))
all_df$observation <- 1:nrow(all_df)
all_df$unique_locus_gene_mutation <- paste(all_df$observation,all_df$locus_tag,all_df$type,all_df$position,all_df$gene_product,all_df$mutation,sep="::")
all_df$locus_gene_mutation <- paste(all_df$locus_tag,all_df$type,all_df$position,all_df$gene_product,all_df$mutation,sep="::")
all_df$sample_names <- paste(all_df$title, all_df$ID,all_df$antibiotic,all_df$biofilm.planktonic,all_df$O2,all_df$day,sep="::")

#write raw output containing everything
write.csv(all_df,file="df_all.csv",row.names = F)

#subtract out reference annotations
all_df_noref <- all_df %>% filter(type!="RA")

write.csv(all_df_noref,file="df_all_noref.csv",row.names = F)

#reference mutations
ref_mutations <- all_df %>% filter(type=="RA")
write.csv(ref_mutations,file="ref_mutations.csv",row.names=F)

#just missing coverage
mc_df <- all_df %>% filter(type=="MC") 
write.csv(mc_df,file="mc_df.csv",row.names = F)

#just new junction evidence
jc_df <- all_df %>% filter(type=="JC") 
write.csv(jc_df,file="jc_df.csv",row.names = F)

####
##creating a data frame that is 'casted'


all_df_norefnojc <- all_df_noref %>% filter(type!="JC")
all_df_norefnomcjc <- all_df_norefnojc %>% filter(type!="MC")
library("tidyverse")
all_df_norefnomcjc <- all_df_norefnomcjc[,c("sample_names","locus_gene_mutation","frequency")] %>% pivot_longer(cols="frequency")

df_snp_t <- all_df_norefnomcjc %>% pivot_wider(names_from=c("sample_names"),values_from="value", values_fill=0)

df_snp_t <- df_snp_t %>% column_to_rownames(var="locus_gene_mutation")

df_snp_t <- df_snp_t[,-1]

#add two columns, count and sum, summarizing number and frequency of each mutation
df_snp_t_sum <- transform(df_snp_t, count=rowSums(df_snp_t!=0.0), sum=rowSums(df_snp_t))

#add column snp_type so we can filter out nonsynonymous
all_mut_type <- all_df %>% select(c(locus_gene_mutation, snp_type,locus_tag)) %>% distinct()

df_snp_t_sum$locus_gene_mutation <- rownames(df_snp_t_sum)

df_snp_t_sum <- left_join(df_snp_t_sum, all_mut_type,by="locus_gene_mutation")
df_snp_t_sum <- df_snp_t_sum %>% column_to_rownames(var="locus_gene_mutation")


#write casted/pivot_wider data frame that can be sorted by mutation frequency or count
write.csv(df_snp_t_sum,file="all_df_noref_t_count.csv")
