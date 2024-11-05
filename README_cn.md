# MigSeq2: 系统发育分析流程

MigSeq2 是一个用于处理和分析高通量测序数据的流程，主要用于系统发育分析。该流程包括从原始FASTQ文件的质量控制到最终的系统发育树构建。

## 功能

- FASTQ文件质量控制和预处理
- FASTQ到FASTA的转换
- 序列重命名
- 使用Stacks进行denovo分析
- PHYLIP文件格式修改
- 使用IQ-TREE进行系统发育树构建

## 依赖

- Python 3.x
- fastp
- Stacks
- IQ-TREE
- pandas
- pysam
- Bio (Biopython)
- bwa

## 安装

1. 克隆此仓库：
   ```bash
   git clone https://github.com/your_username/MigSeq2.git
   cd MigSeq2
   ```

2. 使用提供的Dockerfile构建Docker镜像：
   ```bash
   docker build -t migseq2 -f Dockerfile_stacks .
   ```

3. 启动Docker镜像：
   ```bash
   docker run -it -v /path/to/your/data:/home/app migseq2 /bin/bash
   ```

## 使用方法

### 使用Migseq_2_denovo.py进行de novo分析

基本用法：
```bash
python Migseq_2_denovo.py -i <input_directory> -o <output_directory> [options]

# 示例
python Migseq_2_denovo.py -i raw -o output
```

### 使用Docker运行de novo分析流程

如果您希望使用Docker运行de novo分析流程，可以使用以下命令：

```bash
docker run --rm -v /path/to/your/data:/data migseq2 python Migseq_2_denovo.py -i /data/input_directory -o /data/output_directory [options]
```

在这里，`/path/to/your/data` 是您本地数据的路径，`/data/input_directory` 和 `/data/output_directory` 是容器内的输入和输出目录。

### 使用Migseq_2_mapping.py进行基因组映射

基本用法：
```bash
python Migseq_2_mapping.py --indir <input_directory> --ref_genome <reference_genome> --outdir <output_directory> --fada <forward_adapter> --rada <reverse_adapter> [options]
```

示例：
```bash
python Migseq_2_mapping.py --indir rawdata --ref_genome ./refgenome/RedCoralContigSpades.fasta --outdir output --fada AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --rada AGATCGGAAGAGCGTCGTGTAGGGAAAGAC
```

### 参数说明

| 参数 | 描述 | 默认值 |
|------|------|--------|
| `-i, --indir` | 输入目录，包含原始FASTQ文件 | 必需 |
| `-o, --outdir` | 输出目录 | 必需 |
| `-q, --quality` | 质量阈值 | 20 |
| `-F, --F_remove` | 从正向读取中移除的碱基数 | 0 |
| `-R, --R_remove` | 从反向读取中移除的碱基数 | 0 |
| `-t, --threads` | 线程数 | 1 |
| `-l, --min_len` | 最小长度阈值 | 30 |
| `-m, --min_cov` | 最小覆盖度 | 3 |
| `-M, --max_dist` | 最大距离 | 2 |
| `-N, --max_dist_secondary` | 次要读取的最大距离 | 4 |
| `--r` | r值 | 0.8 |
| `--pop_opts` | 种群选项 | '' |
| `--from_fa` | 从FASTA文件开始（可选） | None |
| `--popmap` | 种群映射文件（可选） | None |

## 输出

脚本将在指定的输出目录中生成以下文件和目录：

1. `processed_fastq/`: 包含处理后的FASTQ文件
2. `fa/`: 包含转换后的FASTA文件
3. `renamed/`: 包含重命名后的FASTA文件
4. `pl/`: 包含Stacks的输出文件
5. `iq/`: 包含IQ-TREE的输出文件，包括系统发育树

## 注意事项

- 确保输入的FASTQ文件命名格式正确（例如：`sample_L001_R1_001.fastq.gz`）
- 根据您的数据特性调整参数，特别是质量阈值、最小长度和覆盖度
- 对于大型数据集，请确保有足够的计算资源和存储空间

## 许可

[在此添加许可信息]

## 联系

[在此添加联系信息]