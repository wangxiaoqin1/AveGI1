#!/bin/bash
# 两阶段基因岛提取流程
# 阶段1: 用int基因筛选阳性菌株 (identity ≥ 90%)
# 阶段2: 对阳性菌株的完成图，用int+fhub提取基因岛 → Prokka → ResFinder

set -e

# ========== 参数 ==========
INT_FASTA="$1"          # int基因序列文件
FHUB_FASTA="$2"         # fhub基因序列文件
GENOMES_DIR="$3"        # 存放所有基因组完成图的目录
OUT_DIR="$4"            # 输出目录
THREADS="${5:-1}"       # 线程数，默认1
IDENTITY="${6:-90}"     # BLAST一致性，默认90
EXPAND="${7:-100}"      # 提取扩展碱基数，默认100

if [ $# -lt 4 ]; then
    echo "用法: $0 <int.fa> <fhub.fa> <genomes_dir> <out_dir> [threads] [identity] [expand]"
    exit 1
fi

mkdir -p "$OUT_DIR"/{positive_strains,islands,prokka,resfinder}

# ========== 阶段1: 筛选阳性菌株 ==========
echo "===== 阶段1: 基于int基因筛选阳性菌株 (identity ≥ ${IDENTITY}%) ====="
positive_list="$OUT_DIR/positive_strains/positive_genomes.list"
> "$positive_list"

for genome in "$GENOMES_DIR"/*.{fa,fasta,fna,fsa} 2>/dev/null; do
    [ -f "$genome" ] || continue
    base=$(basename "$genome" | sed 's/\.[^.]*$//')
    
    # 快速BLAST int基因
    if blastn -query "$INT_FASTA" -db "$genome" -perc_identity "$IDENTITY" -evalue 1e-5 \
              -outfmt 6 | head -1 | grep -q .; then
        echo "  [+] 阳性: $base"
        echo "$genome" >> "$positive_list"
    else
        echo "  [-] 阴性: $base"
    fi
done

pos_count=$(wc -l < "$positive_list")
echo "阶段1完成: 共筛选出 $pos_count 株阳性菌株"
if [ $pos_count -eq 0 ]; then
    echo "无阳性菌株，流程终止"
    exit 0
fi

# ========== 阶段2: 处理阳性菌株 ==========
echo "===== 阶段2: 提取基因岛并注释 ====="
while IFS= read -r genome; do
    base=$(basename "$genome" | sed 's/\.[^.]*$//')
    echo ">>> 处理阳性菌株: $base"

    # 确保是单contig完成图（可选检查）
    if [ $(grep -c "^>" "$genome") -ne 1 ]; then
        echo "    警告: 不是单contig，跳过"; continue
    fi

    # 建立索引
    makeblastdb -in "$genome" -dbtype nucl -parse_seqids >/dev/null 2>&1

    # 定位int基因（取最佳hit）
    int_hit=$(blastn -query "$INT_FASTA" -db "$genome" -perc_identity "$IDENTITY" \
              -evalue 1e-5 -outfmt "6 sseqid sstart send" | head -1)
    [ -z "$int_hit" ] && { echo "    未找到int基因（异常）"; continue; }
    int_contig=$(echo "$int_hit" | cut -f1)
    int_start=$(echo "$int_hit" | cut -f2)
    int_end=$(echo "$int_hit" | cut -f3)
    [ $int_start -gt $int_end ] && { tmp=$int_start; int_start=$int_end; int_end=$tmp; }

    # 定位fhub基因
    fhub_hit=$(blastn -query "$FHUB_FASTA" -db "$genome" -perc_identity "$IDENTITY" \
               -evalue 1e-5 -outfmt "6 sseqid sstart send" | head -1)
    [ -z "$fhub_hit" ] && { echo "    未找到fhub基因"; continue; }
    fhub_contig=$(echo "$fhub_hit" | cut -f1)
    fhub_start=$(echo "$fhub_hit" | cut -f2)
    fhub_end=$(echo "$fhub_hit" | cut -f3)
    [ $fhub_start -gt $fhub_end ] && { tmp=$fhub_start; fhub_start=$fhub_end; fhub_end=$tmp; }

    # 检查同contig
    if [ "$int_contig" != "$fhub_contig" ]; then
        echo "    int和fhub不在同一contig，跳过"; continue
    fi

    # 提取区间（int和fhub之间的区域）
    start=$(( (int_start < fhub_start ? int_start : fhub_start) - EXPAND ))
    end=$(( (int_end > fhub_end ? int_end : fhub_end) + EXPAND ))
    [ $start -lt 1 ] && start=1

    echo "    提取区域: ${int_contig}:${start}-${end}"
    samtools faidx "$genome" "${int_contig}:${start}-${end}" > "$OUT_DIR/islands/${base}_island.fasta" 2>/dev/null
    [ ! -s "$OUT_DIR/islands/${base}_island.fasta" ] && { echo "    提取失败"; continue; }

    # Prokka注释
    prokka "$OUT_DIR/islands/${base}_island.fasta" --outdir "$OUT_DIR/prokka/${base}" \
           --prefix "$base" --cpus "$THREADS" --force >/dev/null 2>&1

    # ResFinder耐药基因
    run_resfinder.py -ifa "$OUT_DIR/islands/${base}_island.fasta" -o "$OUT_DIR/resfinder/${base}" \
                     -t "$THREADS" >/dev/null 2>&1
    if [ -f "$OUT_DIR/resfinder/${base}/ResFinder_results_tab.txt" ]; then
        echo "    耐药基因:"
        cat "$OUT_DIR/resfinder/${base}/ResFinder_results_tab.txt" | sed 's/^/      /'
    fi
done < "$positive_list"

echo "===== 全部完成 ====="
echo "阳性菌株列表: $positive_list"
echo "基因岛序列: $OUT_DIR/islands/"
echo "Prokka结果: $OUT_DIR/prokka/"
echo "ResFinder结果: $OUT_DIR/resfinder/"
