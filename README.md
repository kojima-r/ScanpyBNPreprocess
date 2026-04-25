# ScanpyBNPreprocess

[Tabula Muris Senis](https://tabula-muris-senis.ds.czbiohub.org/) のシングルセル RNA-seq データを前処理し、複数のベイジアンネットワーク (BN) 推定バックエンドに流し込むためのパイプラインです。

対応している入力データ:

| mode      | 説明 |
| --------- | ---- |
| `bbknn`   | FACS と droplet を統合し BBKNN (Batch Balanced k-Nearest Neighbors) でバッチ補正した統合データ |
| `facs`    | FACS (Fluorescence-Activated Cell Sorting) 方式 |
| `droplet` | ドロップレット方式 |

対応している BN 推定バックエンド (5 種):

| 出口          | 規模    | 入力     | 主な用途                                           |
| ------------- | ------- | -------- | -------------------------------------------------- |
| **ingor**     | 大規模  | 連続値   | ブートストラップ構造学習 + ECv で個別ネットワーク |
| **FastBN**    | 大規模  | 離散値   | C++ 実装の Hill-Climb + Tabu (BIC/BDeu/K2)        |
| **MatDNF**    | 小規模  | 離散値   | DNF ベースの Boolean BN (説明可能 NN)              |
| **ilp\_bn**   | 小規模  | 離散値   | ILP (PuLP/CBC) による厳密最適化                    |
| **pgmpy**     | 小規模  | 離散値   | pgmpy の HillClimb / PC / Hybrid                  |

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
                        h5ad
                          │
                          ▼
                  01preprocess.py
                          │
              per-tissue (cells×genes)  ─── 01data_<mode>/
                          │
        ┌─────────────────┼─────────────────────┐
        ▼                 ▼                     ▼
  02resample.py    02resample.py --batched  02pseudo_bulk.py
   (cells×genes)        (層別)               (各臓器1サンプル)
        │                 │                     │
        ▼                 ▼                     ▼
  02data_bbknn/     02data_bbknn_b/         02data_bbknn2/
        │                 │                     │
        │                 │                     │
        ▼ 03transpose     │                     ▼ 03transpose
  02data_bbknn_t/         │                02data_bbknn2_t/
        │                 │                     │
        ▼ 04merge         ▼ 04merge             ▼ 04merge
  03data_bbknn/all.txt    03data_bbknn_b/  03data_bbknn/all2.txt
   (連続値, genes×samples) all.txt          (連続値, genes×samples)
        │                 │                     │
        │                 ▼ 05disc.py           │
        │           03data_bbknn_b/all_disc.txt │
        │           03data_bbknn_b/all_disc_tri.txt
        │                 │                     │
        │                 ▼ 06prep_disc.py      │
        │           03data_bbknn_b/all_disc{,10,100,1000}.tsv
        │           03data_bbknn_b/tissue/<tissue>.tsv
        │                 │                     │
        ▼                 │                     ▼
   ┌─[ ingor ]─┐    ┌───────────────────────────┴──┐
   │ 05run.sh  │    │   discretized BN backends    │
   │ 06run_ecv │    │ ─ FastBN  (大規模)           │
   └───────────┘    │ ─ MatDNF  (小規模)           │
                    │ ─ ilp_bn  (小規模)           │
                    │ ─ pgmpy   (小規模)           │
                    └──────────────────────────────┘
```

---

## 1. Preprocess — `01preprocess.py`

`data/*.h5ad` から「臓器ごとの (cell × highly variable gene)」行列を書き出します。
1 行目はサンプル ID `tissue|age|batch|cell_id` (FACS / droplet では `batch` 抜き)。

```sh
python 01preprocess.py --mode bbknn
# 出力: 01data_bbknn/<tissue>.txt
```

参考ノートブック: `preprocess_facs.ipynb`, `preprocess_bbknn.ipynb` (UMAP までの確認用)

---

## 2. Resample / Pseudo-bulk

### `02resample.py`
ブートストラップサンプリング。各臓器ファイルにつき `N` 個のリサンプル平均を出力します。

```sh
# 通常版
python 02resample.py --input-glob "01data_bbknn/*.txt" -n 10
# 出力: 02data_bbknn/<tissue>.txt

# age|batch ごとに層別 (batched) ※ 離散化系バックエンド向けに必要
python 02resample.py --input-glob "01data_bbknn/*.txt" -n 10 --batched
# 出力: 02data_bbknn_b/<tissue>.txt
```

オプション: `--size {same,root}` `--workers N`

### `02pseudo_bulk.py`
各臓器の全細胞を 1 サンプルに平均化:

```sh
python 02pseudo_bulk.py --input-glob "01data_bbknn/*.txt"
# 出力: 02data_bbknn2/<tissue>.txt
```

---

## 3. Transpose — `03transpose.py`

(cells × genes) を (genes × cells) に転置します。

```sh
python 03transpose.py
# 02data_bbknn/*  → 02data_bbknn_t/
# 02data_bbknn2/* → 02data_bbknn2_t/

python 03transpose.py --input-glob "02data_bbknn_b/*.txt" --out-dir 02data_bbknn_b_t/
```

---

## 4. Merge — `04merge.py`

各臓器ファイルを結合します。

```sh
# 連続値・列方向結合 (ingor 用):
python 04merge.py
#  → 03data_bbknn/all.txt
#  → 03data_bbknn/all2.txt

# 行方向結合 + tissue 列付加 (離散化系バックエンド用):
python 04merge.py --batched
#  → 03data_bbknn_b/all.txt
```

---

## 5. 離散化 — `05disc.py`

`04merge.py --batched` の出力を、各遺伝子ごとに 0.1% / 75% 分位点で binary / ternary に離散化:

```sh
python 05disc.py
# 出力: 03data_bbknn_b/all_disc.txt      (binary,  @name + tissue 付)
#       03data_bbknn_b/all_disc_tri.txt  (ternary, @name + tissue 付)
```

---

## 6. 離散化系バックエンド向けの整形 — `06prep_disc.py`

`05disc.py` の出力から `@name` と `tissue` を落とし、サイズ別の TSV と臓器別 TSV を切り出します。FastBN / MatDNF / ilp\_bn / pgmpy の 4 つの出口がここで生成されたファイルを共通入力として使います。

```sh
python 06prep_disc.py
# 出力 (03data_bbknn_b/ 配下):
#   all_disc.tsv                 全遺伝子・binary
#   all_disc{10,100,1000}.tsv    先頭 N 列のサブセット
#   all_disc_tri{,10,100,1000}.tsv  ternary 版
#   tissue/<tissue>.tsv          臓器別
#   tissue_tri/<tissue>.tsv      臓器別 ternary
```

---

## BN 推定バックエンド (5 種の出口)

### A. ingor — 大規模・連続値

`03data_bbknn/all.txt` (転置・連続値) を入力として、`ingor` のブートストラップ構造学習を行います。

```sh
# 03data_bbknn/all.txt に対して N=10 でブートストラップ
sh 05run.sh 03data_bbknn/all.txt
# 出力: bs_all/result.ing.*

# all2.txt を N=100 で実行 (ラッパー)
sh 05run2.sh
# 出力: bs_all2/result.ing.*

# 各臓器ごとに個別に推定
sh 05run_each.sh
# 出力: bs_<tissue>/result.ing.*
```

ブートストラップ後、ECv で全体ネットワークと臓器ごとの個別ネットワークを得ます:

```sh
sh 06run_ecv_all.sh        # bs_all + 02data_bbknn_t  → ecv_all/
sh 06run_ecv_all2.sh       # bs_all2 + 02data_bbknn2_t → ecv_all2/
```

可視化は `BN_ecv.ipynb` を参照。

### B. FastBN — 大規模・離散値

C++ 実装の Hill-Climb + Tabu (BIC/BDeu/K2)。`06prep_disc.py` が生成した `03data_bbknn_b/all_disc{,1000}.tsv` を直接読み込みます。

```sh
cd FastBN
sh compile.sh                      # g++ -O3 -std=c++17 fast_bn.cpp -o fast_bn
sh run.sh                          # → 03data_bbknn_b/all_disc1000.tsv で構造学習
sh pred.sh                         # 学習済み構造に対する再スコア + edge-importance
```

`fast_bn` 自体の詳細オプション (BDeu/BIC/K2, ブートストラップモード等) は `FastBN/README.md` を参照。

### C. MatDNF — 小規模・離散値 (Boolean DNF)

各遺伝子を target、それ以外を feature として MatDNFClassifier を 1 本ずつ学習し、得られた DNF に登場する変数を親候補として書き出します。

```sh
cd MatDNF
sh run.sh                                       # = python learn_bn_matdnf.py ../03data_bbknn_b/all_disc100.tsv
# 出力: ../03data_bbknn_b/all_disc100_matdnf.json  ({"<gene>": ["<parent>", ...]})
```

MatDNF パッケージ自体の使い方や JAX/CuPy 実装の選択は `MatDNF/README.md` を参照。

### D. ilp\_bn — 小規模・離散値 (ILP / PuLP)

PuLP の CBC ソルバで親集合選択 + 位相順制約による厳密最適化。

```sh
cd ilp_bn
sh run.sh                                # = python bn_ilp_pulp.py ../03data_bbknn_b/all_disc100.tsv
# 出力: ../03data_bbknn_b/all_disc100_bn_structure.json
```

`--max-parents`, `--score {bic,bdeu}`, `--time-limit` など詳細オプションは `bn_ilp_pulp.py --help` を参照。

### E. pgmpy — 小規模・離散値 (HillClimb / PC / Hybrid)

```sh
cd pgmpy
sh run.sh                                # = python learn_bn_pgmpy.py ../03data_bbknn_b/all_disc100.tsv (hybrid+k2)
# 出力: bn_result.{json,bif,png}
```

オプション (`--estimator hc|pc|hybrid`, `--score k2|bdeu|bic`, `--max-indegree`, …) は `learn_bn_pgmpy.py --help` を参照。

---

## XX. AD マウスとの比較

```sh
# データダウンロード
sh AD/run.sh
```

前処理と個別ネットワーク計算は `AD_mouse.ipynb` を参照。

---

## ディレクトリ早見表

| ディレクトリ            | 内容                                                    |
| ----------------------- | ------------------------------------------------------- |
| `data/`                 | 入力 h5ad                                               |
| `01data_<mode>/`        | 臓器ごとの cells×genes 行列                            |
| `02data_bbknn/`         | resample 出力 (cells×genes)                            |
| `02data_bbknn_b/`       | resample --batched 出力 (層別、離散化系の元データ)     |
| `02data_bbknn2/`        | pseudo_bulk 出力                                        |
| `02data_bbknn_t/`       | resample 出力の転置 (genes×cells)                      |
| `02data_bbknn2_t/`      | pseudo_bulk 出力の転置                                  |
| `03data_bbknn/`         | 全臓器を結合した行列 (`all.txt`, `all2.txt`) — ingor 用 |
| `03data_bbknn_b/`       | batched 結合・離散化・サブセット — 離散化系 4 出口の共通入力 |
| `bs_<name>/`            | `ingor` のブートストラップ結果                          |
| `ecv_all/`, `ecv_all2/` | 臓器ごとの ECv 個別ネットワーク                         |
| `FastBN/`, `MatDNF/`, `ilp_bn/`, `pgmpy/` | 各バックエンドの実装                  |
