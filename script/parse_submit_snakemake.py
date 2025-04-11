# VERSION=4.0.0
# mamba create -n RNApipline_env fastp=0.23.4 fastqc=0.12.1 openjdk=20.0.0 pandas=2.0.2 parallel-fastq-dump=0.6.7 python=3 sambamba=1.0 samtools=1.17 sra-tools=3.0.5 subread=2.0.3 hisat2=2.2.1
# mamba activate RNApipline_env
# pip install multiqc snakemake
import pandas as pd
import os
import re
import argparse as arg
import shutil
import yaml
if __name__ == "__main__":
    parser = arg.ArgumentParser(description='Submit rna-seq alignment and quantitation jobs')
    parser.add_argument('-s', '--samples', type=str, help='samples fastq list with full path. eg: /workspace/SRR1.fastq.gz', action="store",default=False)
    parser.add_argument('-o', '--outdir', type=str, help='name of the output directory', action="store")
    parser.add_argument('-g', '--gtf',nargs='?', type=str, help='path of GTF file ', action="store")
    parser.add_argument('-x', '--hisat2_index',nargs='?', type=str, help='HISAT2 index ', action="store")
    parser.add_argument('-d', '--dry_run',help='dry run do not submit', default="F",choices=["T","F"])
    parser.add_argument('-c', '--cores',nargs='?' ,type=int, default=4, help='Number of cores used. default:4')
    parser.add_argument('-t', '--type', type=str, help='model of works', action="store",choices=["quantify","variant","multiqc","variant_part"])
    parser.add_argument('-f', '--fastp_file',nargs='?', type=str, help='fastp_file', action="store")
    parser.add_argument('-m', '--model', type=str,nargs='?', help='paired-end:PE, single-end:SE', action="store", default="PE",choices=["PE","SE"])
    parser.add_argument('-v', '--vcf',nargs='?', type=str, help='vcf file', action="store")
    parser.add_argument('-l', '--lsf',nargs='?', type=str, help='wether to sub to lsf cluster',default="F",choices=["T","F"], action="store")
    parser.add_argument('--slurm',nargs='?', type=str, help='wether to sub to Slrum cluster',default="F",choices=["T","F"], action="store")
    parser.add_argument('--slurm_partition',nargs='?', type=str, help='slurm partition',default="cn", action="store")

    args = parser.parse_args()
    if args.cores < 2:
         Exception("The number of cores must be than 2!")
    USER_ID = os.popen("whoami").read().strip()
    
#     # check update
#     version_list = []
#     with open(os.path.realpath(__file__)) as f:
#         version_list.append(f.readlines()[0].rstrip().split("=")[1])
#     new_file = "https://raw.githubusercontent.com/luotao-hust/Bioinfomatics-scripts/main/snakemake_rnaseq/parse_submit_snakemake.py"
#     download_file = os.path.join("/home",USER_ID,".config","snakemake_rnaseq","parse_submit_snakemake.py")
    
#     bk_file = os.path.join("/home",USER_ID,".config","snakemake_rnaseq","parse_submit_snakemake_bk.py")
    
#     if not os.path.exists(os.path.join("/home",USER_ID,".config","snakemake_rnaseq")):   
#         os.mkdir(os.path.join("/home",USER_ID,".config","snakemake_rnaseq"))
#     os.system("wget -q -4 --timeout=6 --tries=2 %s -O $s "%(new_file,download_file))
    
#     if os.path.isfile(download_file):
#         with open(download_file) as f:
#             version_list.append(f.readlines()[0].rstrip().split("=")[1])
#     else:
#         version_list.append("4.0.0")
#     version_list.sort()
#     if version_list[1] != version_list[0]:
#         while True:        
#             print("New version is available! current version is %s!  new version is %s")
#             flag = input("Do you want to update it (y/n)?")
#             flag = flag.upper()
#             if flag in ["Y","N"]:
#                 if flag == "N":
#                     break
#                 else:
    
    output_directory=args.outdir
    samples = pd.read_csv(args.samples, header=None, sep=";")
    samples[0] = samples[0].apply(str)
    current_path = os.path.split(os.path.realpath(__file__))[0]
    config_path = os.path.join(current_path,"../resource","config.yaml")
    f = open(config_path)
    config = yaml.safe_load(f)
    f.close()
    config["sample"] = {}
    config["fastp_para"] = {}
    config["software"] = {}
    config["software"]["gatk"] = os.path.join(current_path,"../software","gatk-4.4.0.0","gatk-package-4.4.0.0-local.jar")
    config["software"]["picard"] = os.path.join(current_path,"../software","picard.jar")
    if args.type == "quantify":
        snakefile = os.path.join(current_path,"../resource","Snakemake_"+args.model+"_quantify")
        check_path = "_counts.txt"
    elif args.type == "variant":
        snakefile = os.path.join(current_path,"../resource","Snakemake_"+args.model+"_variant")
        check_path = "_fil_indel.vcf"
        config["resources"]["fasta"] = os.path.join(os.path.dirname(args.hisat2_index),"genome.fa")
        config["resources"]["vcf"] = args.vcf
    elif args.type == "variant_part":    
        snakefile = os.path.join(current_path,"../resource","Snakemake_"+args.model+"_variant_part")
        check_path = "_fil_indel.vcf"
        config["resources"]["fasta"] = os.path.join(os.path.dirname(args.hisat2_index),"genome.fa")
        config["resources"]["vcf"] = args.vcf
    elif args.type == "multiqc":
        check_path = "../multiqc_report.html"
        snakefile = os.path.join(current_path,"../resource","Snakemake_"+args.model+"_multiqc")

    config["output_directory"] = output_directory
    config["user"] = USER_ID
    config["resources"]["threads"] = args.cores
    config["software"]["hisat2_jun2bed"] = os.path.join(current_path,"../resource","hisat2_jun2bed.py")
    if args.gtf:
        config["resources"]["gtf"] = args.gtf
        if not os.path.isfile(os.path.join(os.path.dirname(args.gtf),"known-splicesite.txt")):
            command = "hisat2_extract_splice_sites.py {gtf_path} > {known_splicesite_path}".format(
                    gtf_path = args.gtf,
                    known_splicesite_path = os.path.join(os.path.dirname(args.gtf),"known-splicesite.txt")
            )
            os.system(command)
        config["resources"]["known_splicesite"] = os.path.join(os.path.dirname(args.gtf),"known-splicesite.txt")
    if args.hisat2_index:
        config["resources"]["hisat2_index"] = args.hisat2_index

    if args.fastp_file:
        fastp_para = pd.read_csv(args.fastp_file, header=None, sep=";")
        for index,row in fastp_para.iterrows():
            config["fastp_para"][row[0].strip()] = row[1]
    else:
        for index,row in samples.iterrows():
            config["fastp_para"][row[0].strip()] = " "

    # wether submit to cluster
    if args.dry_run != "F":
        dry="-np"
        sub2cluster = " "
    elif args.dry_run == "F" and args.lsf == "T":
        dry=""
        sub2cluster = " --profile lsf "
        shutil.copy(os.path.join(current_path,"../resource","lsf.yaml"), output_directory)
    elif args.dry_run == "F" and args.slurm == "T":
        dry=""
        sub2cluster=" --executor slurm --default-resources slurm_partition=" + args.slurm_partition
    else:
        dry=""
        sub2cluster = " "
        
    job_div = 8 if args.type == "quantify" else 11
    rm_flag = {}
    rm_flag_sra = False
    if args.model == "PE":
        if samples.shape[1] == 2:
            for index,row in samples.iterrows():
                config["sample"][row[0].strip()] = row[1].strip()
                rm_flag_sra = True
        elif samples.shape[1] == 3:
            for index,row in samples.iterrows():
                config["sample"][row[0].strip()] = [row[1].strip(),row[2].strip()]
                if not os.path.isfile( os.path.join(output_directory,row[0].strip(),row[0].strip()+"_1.fastq.gz")):
                    rm_flag_sra = True
                    rm_flag[row[0].strip()] = True
                    os.makedirs(os.path.join(output_directory,row[0].strip()), exist_ok=True)
                    os.system("ln -s {source_file} {target_file}".format(source_file = row[1].strip(),target_file = os.path.join(output_directory,row[0].strip(),row[0].strip()+"_1.fastq.gz")))
                if not os.path.isfile( os.path.join(output_directory,row[0].strip(),row[0].strip()+"_2.fastq.gz")):
                    rm_flag_sra = True
                    os.makedirs(os.path.join(output_directory,row[0].strip()), exist_ok=True)
                    os.system("ln -s {source_file} {target_file}".format(source_file = row[2].strip(),target_file = os.path.join(output_directory,row[0].strip(),row[0].strip()+"_2.fastq.gz")))
        else:
            raise Exception("sample file error!")
    else:
        if samples.shape[1] == 2:
            for index,row in samples.iterrows():
                config["sample"][row[0].strip()] = row[1].strip()
                if row[1].strip().endswith(".sra"):
                    rm_flag_sra = True
                if not os.path.isfile( os.path.join(output_directory,row[0].strip(),row[0].strip()+".fastq.gz")) and not row[1].strip().endswith(".sra"):
                    rm_flag[row[0].strip()] = True
                    os.makedirs(os.path.join(output_directory,row[0].strip()), exist_ok=True)
                    os.system("ln -s {source_file} {target_file}".format(source_file = row[1].strip(),
                                                                                                                            target_file = os.path.join(output_directory,row[0].strip(),row[0].strip()+".fastq.gz")))
        else:
            raise Exception("sample file error!")

    # write the new config yaml file to output directory this is for reproducibility
    new_yaml = output_directory + "/config.yaml"
    with open(new_yaml, 'w') as outfile:
        yaml.dump(config, outfile, default_flow_style=False)
    shutil.copy(snakefile, output_directory)
    shutil.copy(os.path.join(current_path,"../resource","env.yaml"), output_directory)
    shutil.copy(os.path.join(current_path,"../resource","lsf.yaml"), output_directory)
    shutil.copy(os.path.join(current_path,"../resource","env_hisat2.yaml"), output_directory)
    shutil.copy(os.path.join(current_path,"../resource","hisat2_jun2bed.py"), output_directory)
    shutil.copy(os.path.join(current_path,"../resource","count_to_tpm.py"), output_directory)
    
    run_commd = "cd {output_directory} && snakemake {sub2cluster} -s {snakefile} --resources mem_mb=50000 --cores {cores} --jobs {jobs} --configfile={config} -k {dry} ".format(
            output_directory = output_directory,
            sub2cluster = sub2cluster,
            snakefile = os.path.join(output_directory,os.path.split(snakefile)[1]),
            cores = args.cores,
            jobs = args.cores // job_div if args.cores > job_div else 1,
            config = os.path.join(new_yaml),
            dry = dry
    )
    os.system(run_commd)
    for index,row in samples.iterrows():
        if os.path.exists(os.path.join(output_directory,row[0].strip(),row[0].strip()+check_path)):
            print(row[0].strip()+" Done!")
            if row[0].strip() in rm_flag.keys():
                if rm_flag[row[0].strip()] and args.type != "multiqc":
                    os.system("rm -f {fastq_gz} {R1_fastq_gz} {R2_fastq_gz}".format(fastq_gz = os.path.join(output_directory,row[0].strip(),row[0].strip()+".fastq.gz"),
                                                                                    R1_fastq_gz = os.path.join(output_directory,row[0].strip(),row[0].strip()+"_1.fastq.gz"),
                                                                                    R2_fastq_gz = os.path.join(output_directory,row[0].strip(),row[0].strip()+"_2.fastq.gz")))
                if rm_flag_sra:
                    os.system("rm -f {fastq_gz} {R1_fastq_gz} {R2_fastq_gz}".format(fastq_gz = os.path.join(output_directory,row[0].strip(),row[0].strip()+".fastq.gz"),
                                                                                    R1_fastq_gz = os.path.join(output_directory,row[0].strip(),row[0].strip()+"_1.fastq.gz"),
                                                                                    R2_fastq_gz = os.path.join(output_directory,row[0].strip(),row[0].strip()+"_2.fastq.gz")))
        os.system("rm -rf {config} {snakefile} {py_file} {snakemake_env} {env} {env_hisat2} {count_to_tpm} {lsf}".format(
            snakefile = os.path.join(output_directory,os.path.split(snakefile)[1]),
            config = os.path.join(new_yaml),
            py_file = os.path.join(output_directory,"hisat2_jun2bed.py"),
            snakemake_env = os.path.join(output_directory,".snakemake"),
            env = os.path.join(output_directory,"env.yaml"),
            env_hisat2 = os.path.join(output_directory,"env_hisat2.yaml"),
            count_to_tpm = os.path.join(output_directory,"count_to_tpm.py"),
            lsf = os.path.join(output_directory,"lsf.yaml")
        ))
        
