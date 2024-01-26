args <- commandArgs()
file_name = args[6]
running_path = args[7]
system(paste0("mkdir -p ",running_path))
system(paste0("scp luot@10.10.1.247:",paste0("/workspace2/luot/CML/script/gatk_call_mu/",file_name," ",running_path)))
sample_list <- read.csv(paste0(running_path,"/",file_name),sep = ";")
for (i in seq_len(nrow(sample_list))){
    system(paste0("mkdir -p ",running_path,"/",sample_list[i,"Sample_id"]))
    system(paste0("scp luot@10.10.1.247:",sample_list[i,"Path"],"/",sample_list[i,"Sample_id"],"/",sample_list[i,"Sample_id"],".sra ",running_path,"/",sample_list[i,"Sample_id"]))
    # system(paste0("scp luot@10.10.1.247:",sample_list[i,"Path"],"/",sample_list[i,"Sample_id"],"/",sample_list[i,"Sample_id"],".bam.bai ",running_path,"/",sample_list[i,"Sample_id"]))
}
sample_list$bam = file.path(running_path,sample_list$Sample_id,paste0(sample_list$Sample_id,".bam"))
sample_list <- sample_list[,c("Sample_id","bam")]
write.table(sample_list,
            file = paste0(running_path,"/sample_list.txt"),sep = ";",quote = F,row.names = FALSE,col.names = FALSE)
