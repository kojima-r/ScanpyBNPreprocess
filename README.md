# ScanpyBNPreprocess

[Tabula Muris Senis](https://tabula-muris-senis.ds.czbiohub.org/) のシングルセル RNA-seq データを前処理し、`ingor` でベイジアンネットワーク (BN) を推定し、ECv で個別ネットワークを得るためのパイプラインです。

対応している入力データ:

| mode      | 説明 |
| --------- | ---- |
| `bbknn`   | FACS と droplet を統合し BBKNN (Batch Balanced k-Nearest Neighbors) でバッチ補正した統合データ |
| `facs`    | FACS (Fluorescence-Activated Cell Sorting) 方式 |
| `droplet` | ドロップレット方式 |

---

## 0. データのダウンロード

```sh
mkdir -p data
wget --content-disposition https://ndownloader.figshare.com/files/23936555 \
     -O data/tabula-muris-senis-bbknn-processed-official-annotations.h5ad
wget --content-disposition https://ndownloader.figshare.com/files/23937842 \
     -O data/tabula-muris-senis-facs-processed-official-annotations.h5ad
```

`00info.py` は figshare の article メタデータを表示する補助スクリプトです。

---

## パイプライン全体像

```
h5ad ── 01preprocess ── per-tissue matrix
                              │
                              ├── 02resample          ── 02data_bbknn/         (cells×genes, bootstrap)
                              ├── 02resample --batched ── 02data_bbknn_b/       (age|batch ごとに層別)
                              └── 02pseudo_bulk       ── 02data_bbknn2/        (各臓器 1 サンプル)
                                                            │
                                                            ▼
                                                       03transpose             (genes×samples)
                                                            │
                                                            ▼
                                                       04merge                 (全臓器を結合)
                                                            │
                                                            ▼
                                                       (05disc, optional)      ← BN 推定向けに離散化
                                                            │
                                                            ▼
                                                       05run.sh                ← ingor で BN 構造推定
                                                            │
                                                            ▼
                                                       06run_ecv_all.sh        ← 個別ネットワーク (ECv)
```

---

## 1. Preprocess — `01preprocess.py`

`data/*.h5ad` から「臓器ごとの (cell × highly variable gene)」行列を書き出します。
1 行目はサンプル ID `tissue|age|batch|cell_id` (FACS / droplet では `batch` 抜き)。

```sh
python 01preprocess.py --mode bbknn
# 入力 : data/tabula-muris-senis-bbknn-processed-official-annotations.h5ad
# 出力 : 01data_bbknn/<tissue>.txt
```

参考ノートブック: `preprocess_facs.ipynb`, `preprocess_bbknn.ipynb` (UMAP までの確認用)

---

## 2. Resample / Pseudo-bulk

### `02resample.py`
ブートストラップサンプリングを行います。各臓器ファイルにつき `N` 個のリサンプル平均を出力します。

```sh
# 通常版
python 02resample.py --input-glob "01data_bbknn/*.txt" -n 10
# 出力: 02data_bbknn/<tissue>.txt

# age|batch ごとの層別 (batched)
python 02resample.py --input-glob "01data_bbknn/*.txt" -n 10 --batched
# 出力: 02data_bbknn_b/<tissue>.txt
```

オプション:

* `--size {same,root}` — リサンプルサイズを元のサンプル数 (`same`, デフォルト) か √N (`root`) にする
* `--workers N` — 並列ワーカー数 (デフォルト 16)

### `02pseudo_bulk.py`
各臓器の全細胞を 1 サンプルに平均化します。

```sh
python 02pseudo_bulk.py --input-glob "01data_bbknn/*.txt"
# 出力: 02data_bbknn2/<tissue>.txt
```

---

## 3. Transpose — `03transpose.py`

(cells × genes) を (genes × cells) に転置します。

```sh
# 引数なしで両方処理 (後方互換):
#   02data_bbknn/*  → 02data_bbknn_t/
#   02data_bbknn2/* → 02data_bbknn2_t/
python 03transpose.py

# 任意の入出力を指定:
python 03transpose.py --input-glob "02data_bbknn_b/*.txt" --out-dir 02data_bbknn_b_t/
```

---

## 4. Merge — `04merge.py`

各臓器のファイルを結合し、全臓器が入った 1 ファイルにします。

```sh
# 引数なしで両方処理 (後方互換):
#   02data_bbknn_t/*  → 03data_bbknn/all.txt
#   02data_bbknn2_t/* → 03data_bbknn/all2.txt
python 04merge.py

# batched 版 (行方向に結合し tissue 列を付加):
python 04merge.py --batched
# 02data_bbknn_b/* → 03data_bbknn_b/all.txt
```

`--input-glob` と `--out` を渡せば任意の組み合わせも可能です。

---

## 5. BN 構造推定 — `05run.sh`

`ingor` を 10 並列のブートストラップで実行します。

```sh
# 03data_bbknn/all.txt に対して N=10 でブートストラップ
sh 05run.sh 03data_bbknn/all.txt
# 出力: bs_all/result.ing.*

# all2.txt に対して N=100 (ラッパー: 05run2.sh)
sh 05run2.sh
# 出力: bs_all2/result.ing.*

# 各臓器ごとに個別に推定
sh 05run_each.sh
# 出力: bs_<tissue>/result.ing.*
```

### 離散化版 — `05disc.py`
`04merge.py --batched` の出力 (`03data_bbknn_b/all.txt`) を、各遺伝子ごとに 0.1% / 75% 分位点で binary / ternary に離散化します。

```sh
python 05disc.py
# 出力: 03data_bbknn_b/all_disc.txt      (binary)
#       03data_bbknn_b/all_disc_tri.txt  (ternary)
```

---

## 6. 個別ネットワーク (ECv) — `06run_ecv_all.sh`

ブートストラップ結果を集約し全体ネットワークを作ったあと、臓器ごとに ECv で個別ネットワークを推定します。

```sh
# bs_all + 02data_bbknn_t  → ecv_all/  (デフォルト)
sh 06run_ecv_all.sh

# bs_all2 + 02data_bbknn2_t → ecv_all2/ (ラッパー: 06run_ecv_all2.sh)
sh 06run_ecv_all2.sh
```

可視化は `BN_ecv.ipynb` を参照。

---

## XX. AD マウスとの比較

```sh
# データダウンロード
sh AD/run.sh
```

前処理と個別ネットワーク計算は `AD_mouse.ipynb` を参照。

---

## ディレクトリ早見表

| ディレクトリ           | 内容                                                |
| --------------------- | --------------------------------------------------- |
| `data/`               | 入力 h5ad                                           |
| `01data_<mode>/`      | 臓器ごとの cells×genes 行列                         |
| `02data_bbknn/`       | resample 出力 (cells×genes)                         |
| `02data_bbknn_b/`     | resample --batched 出力                             |
| `02data_bbknn2/`      | pseudo_bulk 出力                                    |
| `02data_bbknn_t/`     | resample 出力の転置 (genes×cells)                   |
| `02data_bbknn2_t/`    | pseudo_bulk 出力の転置                              |
| `03data_bbknn/`       | 全臓器を結合した行列 (`all.txt`, `all2.txt`)        |
| `03data_bbknn_b/`     | batched 版の結合行列と離散化版                      |
| `bs_<name>/`          | `ingor` のブートストラップ結果                      |
| `ecv_all/`, `ecv_all2/` | 臓器ごとの ECv 個別ネットワーク                   |
