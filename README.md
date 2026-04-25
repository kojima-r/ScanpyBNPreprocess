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
```

**入力例:** `data/tabula-muris-senis-bbknn-processed-official-annotations.h5ad` (AnnData, ~23 万 cells × ~22k genes)

**出力例:** `01data_bbknn/Aorta.txt` (TSV, cells × highly variable genes)

```
@name                                                  Sox17     Mrpl15     Msc        Rdh10      ...
Aorta|18m|1|A15_B002505_B008538_S15.mm10-plus-5-0-1    0.0       1.110112   0.0        0.141129   ...
Aorta|18m|1|A16_B002505_B008538_S16.mm10-plus-5-0-1    0.0       0.892341   0.0        0.0        ...
...
```

参考ノートブック: `preprocess_facs.ipynb`, `preprocess_bbknn.ipynb` (UMAP までの確認用)

---

## 2. Resample / Pseudo-bulk

**共通の入力例:** `01data_bbknn/Aorta.txt` (上記 §1 の出力)

### `02resample.py`
ブートストラップサンプリング。各臓器ファイルにつき `N` 個のリサンプル平均を出力します。

```sh
# 通常版
python 02resample.py --input-glob "01data_bbknn/*.txt" -n 10
```

**出力例 (通常版):** `02data_bbknn/Aorta.txt`

```
@name   Sox17       Mrpl15      Msc         Rdh10       ...
s0      0.831252    0.408112    0.263675    0.209179    ...
s1      0.827084    0.424911    0.271551    0.174305    ...
...                                                       (-n 行 = 10 行)
```

```sh
# age|batch ごとに層別 (batched) ※ 離散化系バックエンド向けに必要
python 02resample.py --input-glob "01data_bbknn/*.txt" -n 10 --batched
```

**出力例 (--batched):** `02data_bbknn_b/Aorta.txt` (`s<i>|<age>|<batch>` ごとに 1 行)

```
@name        Sox17       Mrpl15      Msc         Rdh10       ...
s0|18m|0     0.0         0.919072    0.0         0.097879    ...
s0|18m|1     0.002299    0.597581    0.0         0.053727    ...
s0|24m|1     1.533844    0.408933    0.001413    0.051832    ...
...
```

オプション: `--size {same,root}` `--workers N`

### `02pseudo_bulk.py`
各臓器の全細胞を 1 サンプルに平均化:

```sh
python 02pseudo_bulk.py --input-glob "01data_bbknn/*.txt"
```

**出力例:** `02data_bbknn2/Aorta.txt` (1 行のみ)

```
@name        Sox17       Mrpl15      Msc         Rdh10       ...
sAorta.txt   0.834489    0.411865    0.267889    0.208933    ...
```

---

## 3. Transpose — `03transpose.py`

(cells × genes) を (genes × cells) に転置します。

```sh
python 03transpose.py
# 02data_bbknn/*  → 02data_bbknn_t/
# 02data_bbknn2/* → 02data_bbknn2_t/
```

**入力例:** `02data_bbknn/Aorta.txt` (上記 §2 の出力, `samples × genes`)
**出力例:** `02data_bbknn_t/Aorta.txt` (`genes × samples`)

```
@name      s0          s1          s2          ...   s9
Sox17      0.831252    0.827084    0.872701    ...   0.836877
Mrpl15     0.408112    0.424911    0.483574    ...   0.373387
Msc        0.263675    0.271551    0.234812    ...   ...
...
```

任意の組み合わせも指定可:

```sh
python 03transpose.py --input-glob "02data_bbknn_b/*.txt" --out-dir 02data_bbknn_b_t/
```

---

## 4. Merge — `04merge.py`

各臓器ファイルを結合します。

### 連続値・列方向結合 (ingor 用)

```sh
python 04merge.py
```

**入力例:** `02data_bbknn_t/Aorta.txt`, `02data_bbknn_t/BAT.txt`, ... (§3 の出力, `genes × samples`)
**出力例:** `03data_bbknn/all.txt` (列名は `<tissue>.txt_<sample>`)

```
@name      Skin.txt_s0    Skin.txt_s1    ...    Aorta.txt_s0    Aorta.txt_s1    ...
Sox17      0.002315       0.000415       ...    0.831252        0.827084        ...
Mrpl15     0.658357       0.660373       ...    0.408112        0.424911        ...
...
```

`02data_bbknn2_t/*` も同時に処理され `03data_bbknn/all2.txt` に書き出されます。

### 行方向結合 + tissue 列付加 (離散化系バックエンド用)

```sh
python 04merge.py --batched
```

**入力例:** `02data_bbknn_b/Aorta.txt`, `02data_bbknn_b/BAT.txt`, ... (§2 --batched の出力)
**出力例:** `03data_bbknn_b/all.txt` (末尾に `tissue` 列)

```
@name        Sox17       Mrpl15      Msc         ...    tissue
s0|18m|0     0.0         0.919072    0.0         ...    Aorta
s0|18m|1     0.002299    0.597581    0.0         ...    Aorta
s0|24m|1     1.533844    0.408933    0.001413    ...    Aorta
s0|18m|0     ...                                        BAT
...
```

---

## 5. 離散化 — `05disc.py`

`04merge.py --batched` の出力を、各遺伝子ごとに 0.1% / 75% 分位点で binary / ternary に離散化:

```sh
python 05disc.py
```

**入力例:** `03data_bbknn_b/all.txt` (§4 --batched の出力)
**出力例:** `03data_bbknn_b/all_disc.txt` (binary, `@name` と `tissue` 列付)

```
@name        Sox17  Mrpl15  Msc  Rdh10  ...    tissue
s0|01m|0     1      1       1    1      ...    Aorta
s0|01m|0     1      1       0    1      ...    BAT
s0|03m|0     0      1       1    1      ...    Brain_Myeloid
...
```

**出力例:** `03data_bbknn_b/all_disc_tri.txt` (ternary, 値は 0/1/2)

```
@name        Sox17  Mrpl15  Msc  Rdh10  ...    tissue
s0|01m|0     2      2       1    2      ...    Aorta
s0|01m|0     2      1       0    1      ...    BAT
...
```

---

## 6. 離散化系バックエンド向けの整形 — `06prep_disc.py`

`05disc.py` の出力から `@name` と `tissue` を落とし、サイズ別の TSV と臓器別 TSV を切り出します。FastBN / MatDNF / ilp\_bn / pgmpy の 4 つの出口がここで生成されたファイルを共通入力として使います。

```sh
python 06prep_disc.py
```

**入力例:** `03data_bbknn_b/all_disc.txt`, `03data_bbknn_b/all_disc_tri.txt` (§5 の出力)

**出力例 (すべて `03data_bbknn_b/` 配下):**

| ファイル                              | 内容                                |
| ------------------------------------- | ----------------------------------- |
| `all_disc.tsv`                        | 全遺伝子・binary                    |
| `all_disc{10,100,1000}.tsv`           | 先頭 N 列のサブセット               |
| `all_disc_tri.tsv` ほか               | ternary 版 (同様にサイズ別)         |
| `tissue/<tissue>.tsv`                 | 臓器別 (binary)                     |
| `tissue_tri/<tissue>.tsv`             | 臓器別 (ternary)                    |

**サンプル:** `03data_bbknn_b/all_disc100.tsv` (header 行 + 0/1 のみ、`@name`/`tissue` 列なし)

```
Sox17  Mrpl15  Msc  Rdh10  Stau2  Gdap1  Pi15  Il17a  Gsta3  Rims1  ...   (100 列)
1      1       1    1      1      1      1     1      1      1      ...
1      1       0    1      1      1      0     0      1      1      ...
...
```

---

## BN 推定バックエンド (5 種の出口)

### A. ingor — 大規模・連続値

`03data_bbknn/all.txt` (転置・連続値) を入力として、`ingor` のブートストラップ構造学習を行います。

```sh
sh 05run.sh 03data_bbknn/all.txt          # N=10 でブートストラップ
sh 05run2.sh                              # all2.txt を N=100 (ラッパー)
sh 05run_each.sh                          # 各臓器ごとに個別に推定
```

**入力例:** `03data_bbknn/all.txt` (§4 の出力, `genes × samples`, 連続値)
**出力例:** `bs_all/result.ing.000002`, `result.ing.000003`, ... (ブートストラップ毎の JSON)

```json
{
  "version":1,
  "node":[
    {"id":0, ...},
    ...
  ]
}
```

`05run_each.sh` は `bs_<tissue>/result.ing.*` を出力します。

ブートストラップ後、ECv で全体ネットワークと臓器ごとの個別ネットワークを得ます:

```sh
sh 06run_ecv_all.sh        # bs_all + 02data_bbknn_t  → result_all.sgn3 / ecv_all/
sh 06run_ecv_all2.sh       # bs_all2 + 02data_bbknn2_t → result_all2.sgn3 / ecv_all2/
```

**入力例:** `bs_all/result.ing.*` + `02data_bbknn_t/<tissue>.txt`
**出力例 (全体):** `result_all.sgn3`, `result_all.txt`
**出力例 (臓器別):** `ecv_all/Aorta.txt`

```
Parent  Child   id  parent.id  child.id  edgeScore  BS.Prob   BS.Direction  bspline                 gain      ECv:s0     ECv:s1     ...
Ptprb   Sox17   1   235        1         0.0888889  0.0888889 1.0000000     0.0004608:1.4816226:... 0.0001    -0.084271  -0.074180  ...
Lims2   Sox17   2   1026       1         0.7111111  0.7111111 1.0000000     0.0008235:1.8957900:... 0.0023    0.372293   -0.176525  ...
...
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

**入力例:** `03data_bbknn_b/all_disc1000.tsv` (§6 の出力, 0/1 表)
**出力例 (`FastBN/` 配下):**

| ファイル                  | 内容                                                  |
| ------------------------- | ----------------------------------------------------- |
| `init_edges.tsv`          | 学習された DAG のエッジ (`u  v` 形式)                 |
| `init_edges_named.tsv`    | 同上、遺伝子名付き                                    |
| `all_counts.tsv`          | CPT 推定用のノードカウント                            |
| `edge_importance.tsv`     | `pred.sh` 実行後、エッジごとの重要度                  |

`fast_bn` 自体の詳細オプション (BDeu/BIC/K2, ブートストラップモード等) は `FastBN/README.md` を参照。

### C. MatDNF — 小規模・離散値 (Boolean DNF)

各遺伝子を target、それ以外を feature として MatDNFClassifier を 1 本ずつ学習し、得られた DNF に登場する変数を親候補として書き出します。

```sh
cd MatDNF
sh run.sh                # = python learn_bn_matdnf.py ../03data_bbknn_b/all_disc100.tsv
```

**入力例:** `03data_bbknn_b/all_disc100.tsv` (§6 の出力)
**出力例:** `03data_bbknn_b/all_disc100_matdnf.json`

```json
{
  "Sox17":  ["Mrpl15", "Rims1"],
  "Mrpl15": ["Gdap1"],
  "Msc":    [],
  "Rdh10":  ["Sox17", "Stau2"],
  ...
}
```

MatDNF パッケージ自体の使い方や JAX/CuPy 実装の選択は `MatDNF/README.md` を参照。

### D. ilp\_bn — 小規模・離散値 (ILP / PuLP)

PuLP の CBC ソルバで親集合選択 + 位相順制約による厳密最適化。

```sh
cd ilp_bn
sh run.sh                # = python bn_ilp_pulp.py ../03data_bbknn_b/all_disc100.tsv
```

**入力例:** `03data_bbknn_b/all_disc100.tsv` (§6 の出力)
**出力例:** `03data_bbknn_b/all_disc100_bn_structure.json`

```json
{
  "status": "Optimal",
  "objective_value": -25491.6,
  "order": {"Sox17": 0.0, "Mrpl15": 9.0, "Msc": 9.0, ...},
  "selected_parent_sets": {
    "Sox17":  [],
    "Mrpl15": ["Pi15", "Rims1"],
    "Msc":    ["Rims1", "Sox17"],
    ...
  },
  "edges": [["Pi15", "Mrpl15"], ["Rims1", "Mrpl15"], ...],
  "score": "bdeu", "ess": 1.0, "time_limit": 3600
}
```

`--max-parents`, `--score {bic,bdeu}`, `--time-limit` など詳細オプションは `bn_ilp_pulp.py --help` を参照。

### E. pgmpy — 小規模・離散値 (HillClimb / PC / Hybrid)

```sh
cd pgmpy
sh run.sh                # = python learn_bn_pgmpy.py ../03data_bbknn_b/all_disc100.tsv (hybrid+k2)
```

**入力例:** `03data_bbknn_b/all_disc100.tsv` (§6 の出力)
**出力例:** `pgmpy/bn_result_edges.json` (有向エッジのリスト)

```json
[
  ["Gdap1", "Sox17"],
  ["Gdap1", "Rims1"],
  ["Gdap1", "Il17a"],
  ["Sox17", "Pi15"],
  ["Rims1", "Sox17"],
  ...
]
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
