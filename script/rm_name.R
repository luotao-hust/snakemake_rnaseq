args <- commandArgs()

file_name = args[6]
file_path = args[7]
sample_list <- read.csv(file_name,sep = ";",header = F)
sample_list$bam = file.path(file_path,sample_list$V1,paste0(sample_list$V1,".bam"))
sample_list <- sample_list[,c("V1","bam")]
write.table(sample_list,
            file = file_name,sep = ";",quote = F,row.names = FALSE,col.names = FALSE)


