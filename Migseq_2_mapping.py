import argparse
import pandas as pd
import configparser
import logging
import os
import subprocess
from subprocess import PIPE
import sys
import pysam
import tempfile
from concurrent.futures import ThreadPoolExecutor


class SampleCodeClass:
    def __init__(
            self, indir, q, fada, rada,
            F_remove, R_remove, t, min_len,
            ref_genome, outdir, r, pop_opts):
        self.indir = indir
        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)
        self.log = os.path.join(outdir, 'log.txt')
        self.setup_logger()
        self.q = q
        self.F_remove = F_remove
        self.R_remove = R_remove
        self.fada = fada
        self.rada = rada
        self.t = t
        self.min_len = min_len
        self.ref_genome = ref_genome
        self.path_to_bbmap = "/usr/local/src/bbmap"
        self.path_to_bwa = "bwa"
        self.path_to_samtools = "samtools"
        self.path_to_stacks = "/usr/local/src/stacks-2.64"
        self.r = r
        self.path_to_IQtree = "iqtree2"
        self.pop_opts = pop_opts

    def make_sample_ini(self):
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
        print(f"csv file'{self.ini}'created.")

    def setup_logger(self, name=None):
        logger = logging.getLogger(name)
        logger.setLevel(logging.DEBUG)
        # creates a file handler that logs messages above DEBUG level
        fh = logging.FileHandler(self.log)
        fh.setLevel(logging.DEBUG)
        fh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s %(funcName)s :\n%(message)s')
        fh.setFormatter(fh_formatter)
        # creates a file handler that logs messages above INFO level
        sh = logging.StreamHandler()
        sh.setLevel(logging.INFO)
        sh_formatter = logging.Formatter(
            '[%(asctime)s] %(levelname)s %(filename)s : %(message)s',
            '%Y-%m-%d %H:%M:%S')
        sh.setFormatter(sh_formatter)
        # add the handlers to logger
        logger.addHandler(fh)
        logger.addHandler(sh)
        self.logger = logger

    def execute_cmd(self, cmd, shell=True):
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
        cmd = 'cutadapt --help'
        self.execute_cmd(cmd)

    def load_ini(self):
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
            cmd += f'--adapter_sequence {self.fada} --adapter_sequence_r2 {self.rada} '
            cmd += '--correction --detect_adapter_for_pe'

            try:
                self.execute_cmd(cmd)
                self.samples[sample] = [outfq1, outfq2]
                self.logger.info(f"样本 '{sample}' 处理完成。")
            except Exception as e:
                self.logger.error(f"处理样本 '{sample}' 时出错: {e}")

        with ThreadPoolExecutor(max_workers=self.t) as executor:
            executor.map(process_sample, self.samples.keys())

    def genome_index(self):
        if not os.path.exists(f"{self.ref_genome}.bwt"):
            cmd = f"{self.path_to_bwa} index {self.ref_genome} "
            self.execute_cmd(cmd)

    def genome_mapping(self):
        os.makedirs(os.path.join(self.outdir, 'sam'), exist_ok=True)
        out_sams = []
        for sample in self.samples.keys():
            # input files
            rpf = self.samples[sample][0]
            rer = self.samples[sample][1]
            # output file
            out_sam = f'{sample}.sam'
            output = os.path.join(self.outdir, 'sam', out_sam)
            cmd = f"{self.path_to_bwa} mem -t {self.t} {self.ref_genome} {self.samples[sample][0]} {self.samples[sample][1]} | {self.path_to_samtools} sort -@2 > {output}"
            subprocess.run(cmd, shell=True, check=True)
            self.samples[sample] = output

    def sort_Sb(self):
        os.makedirs(os.path.join(self.outdir, 'bam'), exist_ok=True)
        out_bams = []
        for sample in self.samples.keys():
            samfile = self.samples[sample]
            bamfile = os.path.join(self.outdir, 'bam', f'{sample}.bam')
            pysam.sort("-o", bamfile, samfile)
            self.samples[sample] = bamfile

    def exclude_multi_unmapped_reads(self):
        excluded = []
        self.samples_for_analysis = {}
        for sample in self.samples:
            cmd = f"{self.path_to_samtools} view -c -F 260 {self.samples[sample]}"
            try:
                stdout, stderr = self.execute_cmd(cmd)
                read_count = int(stdout)
                if read_count == 0:
                    excluded.append(sample)
                else:
                    self.samples_for_analysis[sample] = self.samples[sample]
            except ValueError:
                excluded.append(sample)
        for excluded_sample in excluded:
            print("excluded sample:", excluded_sample)

    def pop_map_out(self, popmap=None):
        if popmap is not None:
            self.popmap = popmap
        else:
            names = list(self.samples_for_analysis.keys())
            names.sort()
            with open(os.path.join(self.outdir, "popmap.tsv"), "w") as f:
                for name in names:
                    f.write(f"{name}\t{name}\n")
            self.popmap = os.path.join(self.outdir, "popmap.tsv")

    def gstacks(self, from_bam=None):
        bamdir = os.path.join(self.outdir, 'bam')
        if from_bam is not None:
            bamdir=from_bam
        cmd = f"{self.path_to_stacks}/gstacks -I {bamdir} -M {self.popmap} -O {self.outdir} -t {self.t}"
        self.execute_cmd(cmd)

    def populations(self):
        cmd = f"{self.path_to_stacks}/populations -P {self.outdir} -M {self.popmap} -R {self.r}  --max-obs-het 0.99 --min-maf 0.01 --vcf --structure --treemix --hwe --plink --genepop --phylip -t {self.t}"
        if len(self.pop_opts) > 0:
            cmd += f' {self.pop_opts}'
        self.execute_cmd(cmd)

    def md_phylip(self):
        # Create a md file and open it for writing
        infile_path = os.path.join(self.outdir, "populations.fixed.phylip") 
        outfile_path = os.path.join(self.outdir,"pop.phylip")
        with open(infile_path,"r") as infile, open(outfile_path, "w") as outfile:
            for line in infile:
                if not line.strip().startswith("#"):
                    outfile.write(line)
        
    def iqtree(self):
        phylip_file = os.path.join(self.outdir, 'pop.phylip')
        out = os.path.join(self.outdir,'iq')
        cmd = f"{self.path_to_IQtree} -s {phylip_file} -m MFP+ASC -bb 1000 --alrt 1000 -nt AUTO -pre {out}"
        self.execute_cmd(cmd)

    def run(self, popmap=None, from_bam=None):
        if from_bam is not None:
            print('从BAM文件开始')
            if popmap is None:
                print('请输入popmap文件')
                sys.exit(1)
            self.popmap = popmap
            self.gstacks(from_bam=from_bam)
            self.populations()
            self.md_phylip()
            self.iqtree()
            return 0
        
        self.make_sample_ini()
        self.load_ini()
        print(self.samples)
        self.process_fastq()
        self.genome_index()
        self.genome_mapping()
        self.sort_Sb()
        self.exclude_multi_unmapped_reads()
        self.pop_map_out()
        self.gstacks()
        self.populations()
        self.md_phylip()
        self.iqtree()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", type=str, default="rawdata", help="Input directory")
    parser.add_argument("--ref_genome", type=str, default="./refgenome/RedCoralContigSpades.fasta", help="Reference genome")
    parser.add_argument("--q", type=int, default= 30, help="quality_value")
    parser.add_argument("--F_remove", type=int, default= 0, help="forward_remove_bp")
    parser.add_argument("--R_remove", type=int, default= 15, help="reverse_remove_bp")
    parser.add_argument("--outdir", type=str, default="outdir", help="Output directory")
    parser.add_argument("-t", type=int, default=4, help="Number of threads")
    parser.add_argument("--fada", type=str, default="GTCAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC", help="fada value")
    parser.add_argument("--rada", type=str, default="CAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAC", help="rada value")
    parser.add_argument("--min_len", type=int, default=80, help="Minimum length")
    parser.add_argument("-r", type=float, default=0.7, help="r value")
    parser.add_argument("--popmap", type=str, default=None, help="Population map file")
    parser.add_argument("--from_bam", type=str, default=None, help="BAM file directory")
    parser.add_argument(
        "-pop_opts", '--populations_options', dest='popo',
        type=str,
        default='',
        help='populations command additional options')
    args = parser.parse_args()

    SCC = SampleCodeClass(
        indir=args.indir,
        ref_genome=args.ref_genome,
        q=args.q,
        F_remove=args.F_remove,
        R_remove=args.R_remove,
        outdir=args.outdir,
        t=args.t,
        fada=args.fada,
        rada=args.rada,
        min_len=args.min_len,
        r=args.r,
        pop_opts=args.popo
    )
    SCC.run(popmap=args.popmap, from_bam=args.from_bam)


if __name__ == '__main__':
    main()