

1. 首先导入Python环境

   ```shell
   source ~/python3_env.sh 
   ```

2. 使用conda创建一个新的配置环境（**这一步不要运行**）

   ```shell
   cd /Share/home/zhangqf8/sunlei/data/RNA-RNA/script/gutmman/sprite-pipeline-master
   conda config --add channels r
   conda config --add channels bioconda
   conda create -p .snakemake/conda \
   	bedtools=2.29.0 trim-galore=0.6.2 fastqc=0.11.8 \
   	bowtie2=2.3.5 samtools=1.9 pysam=0.15.0.1 \
   	numpy=1.17.2 multiqc r-base=3.6.1 \
   	r-ggplot2=3.1.1 r-gplots=3.0.1.1 \
   	r-readr=1.3.1 r-optparse=1.6.2 snakemake
   ```

3. 导入环境

   ```shell
   cd /Share/home/zhangqf8/sunlei/data/RNA-RNA/script/gutmman/sprite-pipeline-master
   conda activate .snakemake/conda/lipan
   ```

4. 运行流程

   ```shell
   snakemake \
      --snakefile Snakefile \
      --cores 10 \
      --configfile config.yaml \
      --wrapper-prefix file:///Share/home/zhangqf8/sunlei/data/RNA-RNA/script/gutmman/sprite-pipeline-master/ \
      --verbose
   ```

   可以提交到集群上：

   ```bash
   bsub -q Z-ZQF -e error -o log -n 20 \
     "snakemake \
        --snakefile Snakefile \
        --cores 20 \
        --configfile config.yaml \
        --wrapper-prefix file:///Share/home/zhangqf8/sunlei/data/RNA-RNA/script/gutmman/sprite-pipeline-master/ \
        --verbose"
   ```

   注意当前目录下的`wrapper/wrapper.py`是从[这里](https://bitbucket.org/snakemake/snakemake-wrappers/src/master/bio/cutadapt/se/wrapper.py)下载的。

