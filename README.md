# RNA-seq分析流程

**更新内容**

2024-03-01: snakemake >=8.0 时，slurm 的运行方式发生改变。

2024-01-26: 添加了countReadPairs 参数

2023-08-08: 弃用apptainer，因为一些不稳定的因素。

2023-04-25: 使用parallel-fastq-dump拆分sra文件。在featurecount 加入TPM值。

2023-02-24: 增加对单端测序的支持。

## 1. **RNA-Seq**

本次镜像中使用的数据处理流程为：   `fastp → hisat2 → featureCount` 同时支持GATK找突变的流程

## 2. 环境的安装

**2.1 安装anaconda**

关于**Miniforge项目：**

[https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge)

This repository holds a minimal installer for [Conda](https://conda.io/) specific to [conda-forge](https://conda-forge.org/). Miniforge allows you to install the conda package manager with the following features pre-configured:

- [conda-forge](https://conda-forge.org/) set as the default (and only) channel.
    - Packages in the base environment are obtained from the [conda-forge channel](https://anaconda.org/conda-forge).
- Optional support for PyPy in place of standard Python interpreter (aka "CPython").
- Optional support for [Mamba](https://github.com/mamba-org/mamba) in place of Conda.
- An emphasis on supporting various CPU architectures (x86_64, ppc64le, and aarch64 including Apple M1).

It can be compared to the [Miniconda](https://docs.conda.io/en/latest/miniconda.html) installer.

**2.2 Mambaforge安装地址**

```bash
# 2023-02-22 时间的版本，如果有新的版本可以自己替换地址
wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash ./Mambaforge-Linux-x86_64.sh
# 下面的形式也是可以安装
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh

mamba init bash
mamba init zsh
# 安装完成后可使用 mamba init bash   或者mamba init zsh 语句将mamba 写入环境变量 
```

**2.3 更改conda源**

```bash
vim ~/.condarc
# 把之前的通通删除，改为以下内容
ssl_verify: true
channels:
  - bioconda
  - daler
  - conda-forge
  - free
  - defaults
show_channel_urls: true
channel_priority: flexible

```

2.4 创建环境

完成mamba的安装之后，使用mamba创建环境

```bash
# 创建环境
mamba create -n RNApipline_env fastp=0.23.4 fastqc=0.12.1 openjdk=20.0.0 pandas=2.0.2 parallel-fastq-dump=0.6.7 python=3 sambamba=1.0 samtools=1.17 sra-tools=3.0.5 subread=2.0.3 hisat2=2.2.1
# 启动环境，后续的运行必须在这个环境中
mamba activate RNApipline_env
# 继续安装软件
pip install multiqc snakemake snakemake-executor-plugin-slurm

```

## 3. 运行脚本的下载

目前尚处于测试阶段，随时可能更新。请自行检测

项目位于github

```bash
# mamba install git
mkdir RNAseq_snakemake && cd RNAseq_snakemake 
git clone https://github.com/luotao-hust/snakemake_rnaseq.git

cd script

```

## 4. 使用说明

**4.1 创建sample_list.txt 文件**

运行镜像前需要准备  sample_list.txt 文件，用于标注样本id和样本文件的地址，同时支持从sra文件或者fastq.gz文件开始。

从fastq.gz文件开始定量时：

4.1.1 双端测序时：（单端测序则没有第三列）

文件中有三列，用分号“;” 分割， 第一列为 SRR样本id  第二列为 SRRid_1.fastq.gz文件的绝对路径  第三列为 SRRid_2.fastq.gz 的绝对路径。！！！！！！！文件的名字的形式必须是 SRRid_1.fastq.gz  和  SRRid_2.fastq.gz  ！！！！！！！！！！！！！！ 

[https://zhuanlan.zhihu.com/p/119472868](https://zhuanlan.zhihu.com/p/119472868) (<-**绝对路径与相对路径**)

如下：

sample_list.txt

`SRR0001;/workspace/luot/SRR0001_1.fastq.gz;/workspace/luot/SRR0001_2.fastq.gz`

`SRR0002;/workspace/luot/SRR0002_1.fastq.gz;/workspace/luot/SRR0002_2.fastq.gz`

4.1.2 从.sra文件开始定量时：

文件中有两列，用分号“;” 分割， 第一列为 SRR样本id  第二列为 .sra 文件的绝对路径

如下：

sample_list.txt

`SRR0001;/workspace/luot/SRR0001.sra`

`SRR0002;/workspace/luot/SRR0002.sra`

**4.2 运行脚本起始运行文件**

./script/parse_submit_snakemake.py

**4.3 脚本接受参数说明**

脚本接受以下参数

`-s      sample_list.txt文件的绝对路径`

`-o      输出文件夹`

`-g      GTF文件路径`

`-x      hisat2 index的路径（注意写法 ）`

[Manual | HISAT2](https://daehwankimlab.github.io/hisat2/manual/)

`-c      使用的cpu核心数，默认使用4个` 非必要参数

`-t      运行的类型，  multiqc;quantify;variant ----其中multiqc;仅仅查看测序质量 ；quantify  定量，常规流程。variant  跑gatk 找突变的流程。`

`-m      输入PE 或者SE 指定是双端测序还是单端测序， PE  双端    SE  单端 默认为PE` 非必要参数

`-f      自定义的fastp质控参数文件 sample_fastp.txt` 非必要参数

`-d      测试模式，填写该参数  -d T 后，将测试脚本能否运行，并输出将运行的代码。`非必要参数

`-v      用于GATK 流程的vcf文件地址`  非必要参数

`-lsf      输入 T 或者 F 是否在lsf集群中运行，默认 F` 非必要参数

`-slrum  输入T或者F 是否在Slurm集群中运行，默认F`  非必要参数

## 5. 运行任务

**5.1 以fastp默认质控参数运行**

需要提供：gtf 文件和HISAT2的索引文件和sample_list.txt 文件

示例代码：

```bash
python parse_submit_snakemake.py -o /data/luot/PRJNA000001 -s /data/luot/sample_file.txt  -g /data/luot/data/RNA_seq_pipline/hisat2_index/Homo_sapiens.GRCh38.106.gtf -x /data/luot/data/RNA_seq_pipline/hisat2_index/genome_snp_tran -c 16 -t quantify -m PE -l T

## 提交后台
nohup python parse_submit_snakemake.py -o /data/luot/PRJNA000001 -s /data/luot/sample_file.txt  -g /data/luot/data/RNA_seq_pipline/hisat2_index/Homo_sapiens.GRCh38.106.gtf -x /data/luot/data/RNA_seq_pipline/hisat2_index/genome_snp_tran -c 16 -t quantify -m PE -l T > rnaseq_snakemake_v0.1.0.quantify.log 2>&1 &!

```

**5.2 自定义fastp质控参数运行**

需要提供：gtf 文件和HISAT2的索引文件和sample_list.txt 文件和自定义的fastp质控参数文件sample_fastp.txt

sample_fastp.txt包括两列，用分号“;” 分割，第一列为样本id，第二列为该id的质控参数。如下：    

`SRR00001;-W 4 -M 30 -r -f 10 -F 10  -l 35`

`SRR00002;-W 4 -M 30 -r -f 10 -F 10  -l 35`

流程分两步，第一，计算出 multiqc的结果，第二，根据multiqc的结果做出质控参数

如果质控参数都是相同的，可用以下代码快速生成，， 替换fastp_paras为对应的参数，比如 `-W 4 -M 30 -r -f 10 -F 10  -l 35`  sample_list.txt 为之前提供的样本信息文件

```bash
# awk -F ";" '{printf("%s;%s\n",$1,"fastp_paras") }' sample_list.txt > sample_fastp.txt
	awk -F ";" '{printf("%s;%s\n",$1,"-W 4 -M 30 -r -f 10 -F 10  -l 35") }' sample_list.txt > sample_fastp.txt # 示例
```

5.2.1 计算出 multiqc的结果

multiqc的文件报告在 -o 指定的输出文件夹中

示例代码：

```bash
python parse_submit_snakemake.py -o /data/luot/PRJNA000001 -s /data/luot/sample_file.txt  -c 4 -t multiqc -m PE -l T

# 提交后台
nohup python parse_submit_snakemake.py -o /data/luot/PRJNA000001 -s /data/luot/sample_file.txt -c 4 -t multiqc -m PE -l T > rnaseq_snakemake_v0.1.0.multiqc.log 2>&1 &!

```

5.2.2 做出质控参数并定量

质控参数文件sample_fastp.txt的形式如上面的描述

示例代码

```bash
python parse_submit_snakemake.py -o /data/luot/PRJNA000001 -s /data/luot/sample_file.txt -f /data/luot/sample_fastp.txt -g /data/luot/data/RNA_seq_pipline/hisat2_index/Homo_sapiens.GRCh38.106.gtf -x /data/luot/data/RNA_seq_pipline/hisat2_index/genome_snp_tran -c 16 -t quantify -m PE -l T
```
