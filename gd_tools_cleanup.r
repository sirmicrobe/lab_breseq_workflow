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
all_snp_subinh <- lapply(temp,read.delim) #create a list with each object in list a data frame of mutation table

# headers we want to keep
headers <- c("aa_new_seq","aa_position","aa_ref_seq","codon_new_seq","codon_number","codon_position","codon_ref_seq","start","end","side_1_position","side_2_position","locus_tag","gene_name","gene_position","gene_product","frequency","mutation_category","new_seq","position","seq_id","size","snp_type","type","title")
#take a look at the first data frame in the list
df_snp <- all_snp_subinh[[1]]
df_snp <- df_snp[,headers] #select only headers we want to keep
head(df_snp)
# convert all data frames in the list to one large data frame
all_snp_df <- as.data.frame(data.table::rbindlist(all_snp_subinh,fill=T),colClasses = c("character"))
all_snp_df <- all_snp_df[,headers] #keep headers we want
all_snp_df <- all_snp_df %>% filter(type!="UN") #remove UN at the end
all_snp_df$mutation <- paste(all_snp_df$aa_ref_seq,all_snp_df$aa_position,all_snp_df$aa_new_seq, sep=":") #create mutation column

head(all_snp_df)
#get rid of incompatible characters in file for R
system(paste("iconv -c -f utf-8 -t ascii",args[2], "> SampleKey_noutf.csv",sep=" "))

sample_key <- read.csv("SampleKey_noutf.csv",header=T,check.names=T)


all_snp_df <- (merge(all_snp_df,sample_key,by="title"))

write.csv(all_snp_df,file="df_snp_noref.csv",row.names = F)