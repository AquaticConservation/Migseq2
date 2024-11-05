# MigSeq2: 系統解析フロー

MigSeq2は、高スループットシーケンシングデータを処理および分析するためのフローで、主に系統解析に使用されます。このフローは、原始FASTQファイルの品質管理から最終的な系統樹の構築までを含みます。

## 機能

- FASTQファイルの品質管理と前処理
- FASTQからFASTAへの変換
- 配列のリネーム
- Stacksを使用したdenovo分析
- PHYLIPファイル形式の変更
- IQ-TREEを使用した系統樹の構築

## 依存関係

- Python 3.x
- fastp
- Stacks
- IQ-TREE
- pandas
- pysam
- Bio (Biopython)

## インストール

1. このリポジトリをクローンします：
   ```bash
   git clone https://github.com/your_username/MigSeq2.git
   cd MigSeq2
   ```

2. 提供されたDockerfileを使用してDockerイメージをビルドします：
   ```bash
   docker build -t migseq2 -f Dockerfile_stacks .
   ```

3. Dockerイメージを起動します：
   ```bash
   docker run -it -v /home/xiafei/Migseq2:/home/app migseq2 /bin/bash
   ```

## 使用方法

基本的な使い方：
```bash
python Migseq_2_denovo.py -i <input_directory> -o <output_directory> -f <forward_adapter> -r <reverse_adapter> [options]

python Migseq_2_denovo.py -i raw -o output -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT 
```

### Dockerを使用して分析フローを実行する

Dockerを使用して分析フローを実行したい場合は、次のコマンドを使用できます：

```bash
docker run --rm -v /path/to/your/data:/data migseq2 python Migseq_2_denovo.py -i /data/input_directory -o /data/output_directory -f <forward_adapter> -r <reverse_adapter> [options]
```

ここで、`/path/to/your/data`はローカルデータのパスで、`/data/input_directory`と`/data/output_directory`はコンテナ内の入力および出力ディレクトリです。

### パラメータの説明

| パラメータ | 説明 | デフォルト値 |
|------|------|--------|
| `-i, --indir` | 入力ディレクトリ、原始FASTQファイルを含む | 必須 |
| `-o, --outdir` | 出力ディレクトリ | 必須 |
| `-q, --quality` | 品質閾値 | 20 |
| `-f, --fada` | 正方向アダプター配列 | 必須 |
| `-r, --rada` | 逆方向アダプター配列 | 必須 |
| `-F, --F_remove` | 正方向リードから削除する塩基数 | 0 |
| `-R, --R_remove` | 逆方向リードから削除する塩基数 | 0 |
| `-t, --threads` | スレッド数 | 1 |
| `-l, --min_len` | 最小長さ閾値 | 30 |
| `-m, --min_cov` | 最小カバレッジ | 3 |
| `-M, --max_dist` | 最大距離 | 2 |
| `-N, --max_dist_secondary` | 二次リードの最大距離 | 4 |
| `--r` | r値 | 0.8 |
| `--pop_opts` | 集団オプション | '' |
| `--from_fa` | FASTAファイルから開始（オプション） | None |
| `--popmap` | 集団マッピングファイル（オプション） | None |

## 出力

スクリプトは指定された出力ディレクトリに以下のファイルとディレクトリを生成します：

1. `processed_fastq/`: 処理されたFASTQファイルを含む
2. `fa/`: 変換されたFASTAファイルを含む
3. `renamed/`: リネームされたFASTAファイルを含む
4. `pl/`: Stacksの出力ファイルを含む
5. `iq/`: IQ-TREEの出力ファイルを含む、系統樹を含む

## 注意事項

- 入力のFASTQファイルの命名形式が正しいことを確認してください（例：`sample_L001_R1_001.fastq.gz`）
- データの特性に応じてパラメータを調整してください、特に品質閾値、最小長さ、カバレッジ
- 大規模なデータセットの場合、十分な計算リソースとストレージスペースがあることを確認してください

## ライセンス

[ここにライセンス情報を追加]

## 連絡先

[ここに連絡先情報を追加]
