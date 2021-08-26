#!/usr/bin/env Rscript --vanilla

#name: calculate mapping percentages from breseq directory
#built for output.gd files based on line numbers in breseq v0.35.7

#input: full path to directory containing all breseq subdirectories

library("tidyverse")


#command line inputs
args <- commandArgs(trailingOnly=TRUE)

# Find all .gd files
files <- dir(args[1], recursive=TRUE, full.names=TRUE, pattern="\\output.gd$")

# Make a function to process each file
processFile <- function(files) {
  #read just the lines with mapping info into file
  df <- as.data.frame(read_lines(files[1],skip=9,n_max=4))
  #fix the data frame into two columns
  df2<- separate(df,1,c("label","reads"),"\t")
  #divide mapped by total reads
  as.numeric(df2[3,2])/as.numeric(df2[1,2])
}


# Apply the function to all files.
result <- (sapply(files, processFile))
#write results to tibble
ddf<- tibble(x=as.character(names(result)),y=result)
#filter out redundant values due to "data" subdirectory
ddf <- ddf %>% filter(!str_detect(x,"data"))
#write a csv file
setwd(args[1])
write.csv(ddf,"mapping_percentages.csv",row.names = F)

