# 导入所需的库
import argparse
import pandas as pd
import configparser
import logging
import os
import subprocess
from subprocess import PIPE
import sys
import pysam
import gzip
from Bio import SeqIO
import shutil
from concurrent.futures import ThreadPoolExecutor
import time
import cProfile
import pstats


class SampleCodeClass:
    def __init__(
            self, indir, q,
            F_remove, R_remove, t, min_len, outdir,
            m, M, N, r, pop_opts):
        # 初始化类的属性
        if not os.path.exists(indir):
            raise FileNotFoundError(f"输入目录 '{indir}' 不存在。请检查路径。")
        self.indir = indir  # 输入目录
        self.outdir = outdir  # 输出目录
        os.makedirs(outdir, exist_ok=True)  # 创建输出目录
        self.log = os.path.join(outdir, 'log.txt')  # 日志文件路径
        self.setup_logger()  # 设置日志记录器
        self.q = q  # 质量值阈值
        self.F_remove = F_remove  # 前向读取移除的碱基数
        self.R_remove = R_remove  # 反向读取移除的碱基数
        self.t = t  # 线程数
        self.min_len = min_len  # 最小长度阈值
        self.path_to_bbmap = "/usr/local/src/bbmap"  # BBMap路径
        self.path_to_bwa = "bwa"  # BWA路径
        self.path_to_samtools = "samtools"  # Samtools路径
        self.path_to_bin = "/usr/local/bin/"  # 二进制文件路径
        self.m = m  # 最小覆盖深度
        self.M = M  # 最大距离
        self.N = N  # 次要读取对齐的最大距离
        self.r = r  # r值
        self.pop_opts = pop_opts  # 种群选项
        self.path_to_IQtree= "iqtree2"  # IQ-TREE路径

    def make_sample_ini(self):
        # 创建样本配置文件
        name_list = []
        ini_list = os.listdir(self.indir)
        for name in ini_list:
            if name.endswith("R1_001.fastq.gz"):
                sample_name = name.split('_L001_R1_')[0]
                r1_file = name
                r2_file = sample_name + '_L001_R2_' + name.split('_L001_R1_')[1]
                name_list.append([sample_name, r1_file, r2_file])
        df = pd.DataFrame(name_list, columns=["#Sample_Name", "R1_File", "R2_File"])
        csv = "samples.ini"
        df.to_csv(csv, index=False)
        self.ini = csv
        print(f"csv文件'{self.ini}'已创建。")

    def setup_logger(self, name=None):
        # 设置日志记录器
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        # 创建一个文件处理器，记录DEBUG级别以上的消息
        fh = logging.FileHandler(self.log)
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s %(funcName)s :\n%(message)s')
        fh.setFormatter(fh_formatter)
        # 创建一个流处理器，记录INFO级别以上的消息
        sh = logging.StreamHandler()
        sh.setLevel(logging.INFO)
        sh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s : %(message)s',
            '%Y-%m-%d %H:%M:%S')
        sh.setFormatter(sh_formatter)
        # 将处理器添加到日志记录器
        logger.addHandler(fh)
        logger.addHandler(sh)
        self.logger = logger

    def execute_cmd(self, cmd, shell=True):
        # 执行命令并记录日志
        self.logger.info('[cmd] {}'.format(cmd))
        proc = subprocess.run(
            cmd, shell=shell, stdout=PIPE,
            stderr=PIPE, text=True)
        stdout = proc.stdout
        stderr = proc.stderr
        if len(stdout) > 1:
            self.logger.info(stdout)
        if len(stderr) > 1:
            self.logger.info(stderr)
        return stdout, stderr

    def exec_cutadapt_help(self):
        # 执行cutadapt帮助命令
        cmd = 'cutadapt --help'
        self.execute_cmd(cmd)

    def load_ini(self):
        # 加载配置文件
        with open(self.ini, 'r') as ini:
            lines = ini.readlines()
        self.samples = {}
        for line in lines:
            if line.startswith('#'):
                continue
            data = [i.strip() for i in line.split(sep=',')]
            fq1 = os.path.join(self.indir, data[1])
            fq2 = os.path.join(self.indir, data[2])
            self.samples[data[0]] = [fq1, fq2]

    def process_fastq(self):
        processed = os.path.join(self.outdir, 'processed_fastq')
        os.makedirs(processed, exist_ok=True)

        def process_sample(sample):
            infq1 = self.samples[sample][0]
            infq2 = self.samples[sample][1]

            # 检查输入文件是否存在
            if not os.path.isfile(infq1):
                self.logger.error(f"输入文件 '{infq1}' 不存在。")
                return
            if not os.path.isfile(infq2):
                self.logger.error(f"输入文件 '{infq2}' 不存在。")
                return

            outfq1 = os.path.join(processed, f'{sample}_R1.fastq.gz')
            outfq2 = os.path.join(processed, f'{sample}_R2.fastq.gz')

            cmd = f'fastp -i {infq1} -I {infq2} -o {outfq1} -O {outfq2} '
            cmd += f'-q {self.q} -l {self.min_len} '
            cmd += f'-f {self.F_remove} -F {self.R_remove} '
            cmd += f'-w {self.t} --json {os.path.join(processed, f"{sample}_fastp.json")} '
            cmd += f'--html {os.path.join(processed, f"{sample}_fastp.html")} '
            cmd += '--correction --detect_adapter_for_pe'

            try:
                self.execute_cmd(cmd)
                self.samples[sample] = [outfq1, outfq2]
                self.logger.info(f"样本 '{sample}' 处理完成。")
            except Exception as e:
                self.logger.error(f"处理样本 '{sample}' 时出错: {e}")

        with ThreadPoolExecutor(max_workers=self.t) as executor:
            executor.map(process_sample, self.samples.keys())

    def convert_fastq_to_fasta(self):
        # 将FASTQ文件转换为FASTA格式
        fa = os.path.join(self.outdir, 'fa')
        os.makedirs(fa, exist_ok=True)
        print(f"Directory {fa} created")
        for sample in self.samples.keys():
            fastq_file_r1 = self.samples[sample][0]
            fastq_file_r2 = self.samples[sample][1]
            fasta_file_r1 = os.path.join(fa, f'{sample}.1.fasta')
            fasta_file_r2 = os.path.join(fa, f'{sample}.2.fasta')

            # 使用gzip读取FASTQ文件
            with gzip.open(fastq_file_r1, 'rt') as fq1, gzip.open(fastq_file_r2, 'rt') as fq2:
                SeqIO.convert(fq1, "fastq", fasta_file_r1, "fasta")
                SeqIO.convert(fq2, "fastq", fasta_file_r2, "fasta")
            
            self.samples[sample] = [fasta_file_r1, fasta_file_r2]
        self.inputfasta = [file for sublist in self.samples.values() for file in sublist]

    def rename(self):
        # 重命名FASTA文件中的序列
        # 定义重命名文件的目录
        renamed = os.path.join(self.outdir, 'renamed')
        os.makedirs(renamed, exist_ok=True)

        for sample, file_paths in self.samples.items():
            for i, file_path in enumerate(file_paths):
                # 为重命名的文件创建新的文件名
                renamed_filename = f'{sample}.{i+1}.fasta.gz'
                renamed_file_path = os.path.join(renamed, renamed_filename)
            
                # 打开原始FASTA文件
                with pysam.FastxFile(file_path) as fh_in, gzip.open(renamed_file_path, 'wt') as fh_out:
                    for entry in fh_in:
                        # 替换每个序列的名称
                        entry.name = f'{sample}:{entry.name}/{i+1}'
                        print('>' + entry.name, entry.sequence, sep='\n', file=fh_out)
                    
                # 更新samples字典中的文件路径
                self.samples[sample][i] = renamed_file_path
            
    def pop_map_out(self, popmap=None):
        # 创建或使用种群映射文件
        if popmap is not None:
            self.popmap = popmap
        else:
            names = list(self.samples.keys())
            names.sort()
            with open(os.path.join(self.outdir, "popmap.tsv"), "w") as f:
                for name in names:
                    f.write(f"{name}\t{name}\n")
            self.popmap = os.path.join(self.outdir, "popmap.tsv")
           
    def pl_stacks(self, from_fa=None):
        # 运行Stacks的denovo_map.pl脚本
        renamed = os.path.join(self.outdir, 'renamed')
        pl = os.path.join(self.outdir, 'pl')
        os.makedirs(pl, exist_ok=True)
        if from_fa is not None:
            renamed = from_fa
        cmd = f'denovo_map.pl -M 2 -T {self.t} -o {pl} --popmap {self.popmap} --samples {renamed} --paired -N 4 -X "ustacks:-m 3 --force-diff-len" -X "populations: -M {self.popmap} -R {self.r}  --max-obs-het 0.99 --write-single-snp --min-maf 0.01 --vcf --structure --plink --treemix --phylip -t {self.t} {self.pop_opts}"'
        self.execute_cmd(cmd)
        # ustacks选项说明
        #-f  输入文件路径
        #-i  此样本的唯一整数ID
        #-o  写入结果的输出路径
        #-M  允许的堆栈之间的最大距离（以核苷酸为单位）（默认2）
        #-m  创建堆栈所需的最小覆盖深度（默认3）
        #-N  允许将次要读取对齐到主要堆栈的最大距离（默认：M + 2）
        #-p  启用并行执行，使用num_threads线程
        #-t  输入文件类型。支持的类型：fasta, fastq, gzfasta, 或 gzfastq（默认：猜测）
        #--name  样本的名称（默认：输入文件名减去后缀）
        #-R  保留未使用的读取
        #-H  禁用从次要读取调用单倍型

    def md_phylip(self):
        # 修改PHYLIP文件格式
        # 创建一个md文件并打开它以进行写入
        infile_path = os.path.join(self.outdir,"pl", "populations.fixed.phylip") 
        outfile_path = os.path.join(self.outdir,"pl","pop.phy")
        with open(infile_path,"r") as infile, open(outfile_path, "w") as outfile:
            for line in infile:
                if not line.strip().startswith("#"):
                    outfile.write(line)

    def iqtree(self):
        # 运行IQ-TREE进行系统发育分析
        phy_file = os.path.join(self.outdir, "pl", 'pop.phy')
        outpre = os.path.join(self.outdir, 'iq')
        cmd = f"{self.path_to_IQtree} -s {phy_file} -m MFP+ASC -bb 1000 --alrt 1000 -nt AUTO -pre {outpre}"
        self.execute_cmd(cmd)

    def profile_execution(self, func):
        def wrapper(*args, **kwargs):
            profiler = cProfile.Profile()
            start_time = time.time()
            result = profiler.runcall(func, *args, **kwargs)
            end_time = time.time()
            
            stats = pstats.Stats(profiler)
            stats.sort_stats('cumulative')
            stats.print_stats()
            
            print(f"执行时间: {end_time - start_time:.2f} 秒")
            return result
        return wrapper

    @profile_execution
    def run(self, popmap=None, from_fa=None):
        try:
            if from_fa is not None:
                print('从FASTA文件开始')
                if popmap is None:
                    print('请输入popmap文件')
                    sys.exit(1)
                self.popmap = popmap
                self.pl_stacks(from_fa=from_fa)
                self.md_phylip()
                self.iqtree()
                return 0

            self.make_sample_ini()
            self.load_ini()
            self.process_fastq()  # 这个方法现在包含了质量控制、接头去除和读取修复
            self.convert_fastq_to_fasta()
            self.rename()
            self.pop_map_out()
            self.pl_stacks()
            self.md_phylip()
            self.iqtree()
        except Exception as e:
            self.logger.error(f"运行分析流程时出错: {e}")
            sys.exit(1)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', required=True, help='Input directory')
    parser.add_argument('-o', '--outdir', required=True, help='Output directory')
    parser.add_argument('-q', '--quality', type=int, default=20, help='Quality threshold')
    parser.add_argument('-F', '--F_remove', type=int, default=0, help='Bases to remove from forward reads')
    parser.add_argument('-R', '--R_remove', type=int, default=0, help='Bases to remove from reverse reads')
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('-l', '--min_len', type=int, default=30, help='Minimum length threshold')
    parser.add_argument('-m', '--min_cov', type=int, default=3, help='Minimum coverage')
    parser.add_argument('-M', '--max_dist', type=int, default=2, help='Maximum distance')
    parser.add_argument('-N', '--max_dist_secondary', type=int, default=4, help='Maximum distance for secondary reads')
    parser.add_argument('--r', type=float, default=0.8, help='r value')
    parser.add_argument('--pop_opts', default='', help='Population options')
    parser.add_argument('--from_fa', help='Start from fasta file (optional)')
    parser.add_argument('--popmap', help='Population map file (optional)')
    return parser.parse_args()

def main():
    # 解析命令行参数
    args = get_args()

    # 创建SampleCodeClass实例
    sample_code = SampleCodeClass(
        args.indir,
        args.quality,
        args.F_remove,
        args.R_remove,
        args.threads,
        args.min_len,
        args.outdir,
        args.min_cov,
        args.max_dist,
        args.max_dist_secondary,
        args.r,
        args.pop_opts
    )

    # 运行分析流程
    sample_code.run(popmap=args.popmap, from_fa=args.from_fa)


if __name__ == '__main__':
    main()


