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
              per-tissue (cells×genes)  ─── data01_<mode>/
                          │
        ┌─────────────────┼─────────────────────┐
        ▼                 ▼                     ▼
  02resample.py --level <L>                  02pseudo_bulk.py --level <L>
   <src>=r                                    <src>=p
        │                                          │
        ▼                                          ▼
  data02_bbknn_r_<L>/                       data02_bbknn_p_<L>/
        │                                          │
        ├── 連続値経路 (col-merge) ────────────────┤
        │                                          │
        ▼ 03transpose --source <src> --level <L>   ▼
  data02_bbknn_<src>_<L>_t/                  data02_bbknn_<src>_<L>_t/
        │                                          │
        ▼ 04merge --transposed --source <src> --level <L>
  data03_bbknn_<src>_<L>_t/all.txt           data03_bbknn_<src>_<L>_t/all.txt
                          │
                          ▼
        ┌──────────────────────────────────────┐
        │  continuous BN backends (transposed) │
        └──────────────────────────────────────┘

        └── 離散化経路 (row-merge / disc) ──┐
                                            ▼
                              04merge --source <src> --level <L>
                              data03_bbknn_<src>_<L>/all.txt
                                            │
                                            ▼ 05disc.py
                              data03_bbknn_<src>_<L>/all_disc.txt
                              data03_bbknn_<src>_<L>/all_disc_tri.txt
                                            │
                                            ▼ 06prep_disc.py
                              data04_bbknn_<src>_<L>_disc/all_disc{,10,100,1000}.tsv
                              data04_bbknn_<src>_<L>_disc/tissue/<tissue>.tsv
                                            │
                                            ▼
                              ┌─────────────────────────────┐
                              │   discretized BN backends   │
                              └─────────────────────────────┘

  ※ 05disc.py / 06prep_disc.py は --transposed を渡すと
    data03_bbknn_<src>_<L>_t/all.txt も入力として受け付けます。
    05disc は同ディレクトリに、06prep_disc は
    data04_bbknn_<src>_<L>{,_t}_disc/ に書き出します。
```

---

## 1. Preprocess — `01preprocess.py`

`data/*.h5ad` から「臓器ごとの (cell × highly variable gene)」行列を書き出します。
1 行目はサンプル ID `tissue|age|batch|cell_id` (FACS / droplet では `batch` 抜き)。

```sh
python 01preprocess.py --mode bbknn
```

**入力例:** `data/tabula-muris-senis-bbknn-processed-official-annotations.h5ad` (AnnData, ~23 万 cells × ~22k genes)

**出力例:** `data01_bbknn/Aorta.txt` (TSV, cells × highly variable genes)

```
@name                                                  Sox17     Mrpl15     Msc        Rdh10      ...
Aorta|18m|1|A15_B002505_B008538_S15.mm10-plus-5-0-1    0.0       1.110112   0.0        0.141129   ...
Aorta|18m|1|A16_B002505_B008538_S16.mm10-plus-5-0-1    0.0       0.892341   0.0        0.0        ...
...
```

参考ノートブック: `preprocess_facs.ipynb`, `preprocess_bbknn.ipynb` (UMAP までの確認用)

---

## 2. Resample / Pseudo-bulk

**共通の入力例:** `data01_bbknn/Aorta.txt` (上記 §1 の出力)

### `02resample.py`
ブートストラップサンプリング。`--level` で層別の粒度を指定します (`02pseudo_bulk.py` と同じ `LEVEL_DEPTH = {"tissue": 1, "age": 2, "batch": 3}` 規約 — `@name` を `|` で分割した先頭 N トークンを層別キーに使う)。
出力先の既定値は `data02_bbknn_r_<level>/`、1 列目は `s<i>|<層別キー>` の形。

| `--level` | 層別キー (1 列目)                              | 既定の出力先              |
| --------- | ---------------------------------------------- | ------------------------- |
| `tissue`  | `s<i>\|<tissue>` (各臓器内で 1 プール)         | `data02_bbknn_r_tissue/`  |
| `age`     | `s<i>\|<tissue>\|<age>`                        | `data02_bbknn_r_age/`     |
| `batch`   | `s<i>\|<tissue>\|<age>\|<batch>`               | `data02_bbknn_r_batch/`   |

```sh
python 02resample.py --input-glob "data01_bbknn/*.txt" --level tissue -n 10
python 02resample.py --input-glob "data01_bbknn/*.txt" --level age    -n 10
python 02resample.py --input-glob "data01_bbknn/*.txt" --level batch  -n 10  # 離散化系バックエンド向け
```

**出力例 (`--level tissue`, `-n 2`):** `data02_bbknn_r_tissue/Aorta.txt`

```
@name       Sox17       Mrpl15      Msc         Rdh10       ...
s0|Aorta    0.828619    0.436943    ...
s1|Aorta    0.797671    0.417907    ...
```

**出力例 (`--level age`, `-n 2`):** `data02_bbknn_r_age/Aorta.txt`

```
@name             Sox17       Mrpl15      ...
s0|Aorta|18m      0.507557    0.429221    ...
s0|Aorta|24m      1.347382    0.413791    ...
s0|Aorta|3m       0.731349    0.405789    ...
s1|Aorta|18m      0.460838    0.476399    ...
s1|Aorta|24m      1.573143    0.389519    ...
s1|Aorta|3m       0.745007    0.355852    ...
```

**出力例 (`--level batch`, `-n 2`):** `data02_bbknn_r_batch/Aorta.txt`

```
@name               Sox17       Mrpl15      ...
s0|Aorta|18m|1      0.521136    0.443637    ...
s0|Aorta|24m|1      1.540964    0.389432    ...
s0|Aorta|3m|1       0.667663    0.388408    ...
s1|Aorta|18m|1      0.487469    0.428553    ...
...
```

オプション: `--size {same,root}` `--workers N`

### `02pseudo_bulk.py`
`@name` (`tissue|age|batch|cell_id`) のどの階層で平均化するかを `--level` で選択します。
出力先は `data02_bbknn_p_<level>/`、1 列目は集約キーそのもの。

| `--level` | 集約キー (1 列目)                     | 既定の出力先              |
| --------- | ------------------------------------- | ------------------------- |
| `tissue`  | `<tissue>` (1 行 / ファイル)          | `data02_bbknn_p_tissue/`   |
| `age`     | `<tissue>\|<age>`                     | `data02_bbknn_p_age/`      |
| `batch`   | `<tissue>\|<age>\|<batch>`            | `data02_bbknn_p_batch/`    |

```sh
python 02pseudo_bulk.py --input-glob "data01_bbknn/*.txt" --level tissue
python 02pseudo_bulk.py --input-glob "data01_bbknn/*.txt" --level age
python 02pseudo_bulk.py --input-glob "data01_bbknn/*.txt" --level batch
```

**出力例 (`--level tissue`):** `data02_bbknn_p_tissue/Aorta.txt` (1 行)

```
@name    Sox17       Mrpl15      Msc         Rdh10       ...
Aorta    0.834489    0.411865    0.267889    0.208933    ...
```

**出力例 (`--level age`):** `data02_bbknn_p_age/Aorta.txt`

```
@name        Sox17       Mrpl15      Msc         Rdh10       ...
Aorta|18m    0.508345    0.441201    0.644212    0.248757    ...
Aorta|24m    1.474089    0.428077    0.001429    0.086646    ...
Aorta|3m     0.724628    0.376614    0.106055    0.249393    ...
```

**出力例 (`--level batch`):** `data02_bbknn_p_batch/Aorta.txt`

```
@name          Sox17       Mrpl15      Msc         Rdh10       ...
Aorta|18m|1    0.508345    0.441201    0.644212    0.248757    ...
Aorta|24m|1    1.474089    0.428077    0.001429    0.086646    ...
Aorta|3m|1     0.724628    0.376614    0.106055    0.249393    ...
```

> Note: tissue 階層の値は cell 数で重み付けされた平均で、age/batch 階層の値の単純平均にはなりません。

`--out-dir` を渡せば任意のパスに書き出せます。

---

## 3. Transpose — `03transpose.py`

(cells × genes) を (genes × cells) に転置します。`--source {r,p}` / `--level {tissue,age,batch}` で入出力ディレクトリが自動的に決まります (上書きしたい場合のみ `--input-glob`/`--out-dir`)。

```sh
# resample (--level tissue) を転置
python 03transpose.py --source r --level tissue
#   data02_bbknn_r_tissue/*  →  data02_bbknn_r_tissue_t/

# pseudo_bulk (--level age) を転置
python 03transpose.py --source p --level age
#   data02_bbknn_p_age/*  →  data02_bbknn_p_age_t/
```

**入力例:** `data02_bbknn_r_tissue/Aorta.txt` (上記 §2 の出力, `samples × genes`)
**出力例:** `data02_bbknn_r_tissue_t/Aorta.txt` (`genes × samples`)

```
@name      s0          s1          s2          ...   s9
Sox17      0.875226    0.827113    0.872701    ...   0.836877
Mrpl15     0.433992    0.426842    0.483574    ...   0.373387
Msc        0.310729    0.352021    0.234812    ...   ...
...
```

---

## 4. Merge — `04merge.py`

各臓器ファイルを結合します。

### 行方向結合 + tissue 列付加 (既定 / 離散化系バックエンド用)

既定では `data02_bbknn_<src>_<L>/` (転置前 = `samples × genes`) を行方向に結合し、各サンプルの (tissue, age, batch) を残します。

```sh
python 04merge.py --source r --level batch
#   data02_bbknn_r_batch/*  →  data03_bbknn_r_batch/all.txt
```

**入力例:** `data02_bbknn_r_batch/Aorta.txt`, `data02_bbknn_r_batch/BAT.txt`, ... (§2 --level batch の出力)
**出力例:** `data03_bbknn_r_batch/all.txt` (末尾に `tissue` 列)

```
@name              Sox17       Mrpl15      Msc         ...    tissue
s0|Aorta|18m|1     0.521136    0.443637    0.662757    ...    Aorta
s0|Aorta|24m|1     1.540964    0.389432    0.000353    ...    Aorta
s0|Aorta|3m|1      0.667663    0.388408    0.089176    ...    Aorta
s0|BAT|18m|1       ...                                        BAT
...
```

### 連続値・列方向結合 (`--transposed`, ingor 用)

`--transposed` を渡すと §3 の転置出力 (`data02_bbknn_<src>_<L>_t/`, `genes × samples`) を列方向に結合します。`--source` / `--level` で入出力が自動決定されます。

```sh
# resample (--level tissue) を結合
python 04merge.py --transposed --source r --level tissue
#   data02_bbknn_r_tissue_t/*  →  data03_bbknn_r_tissue_t/all.txt

# pseudo_bulk (--level age) を結合
python 04merge.py --transposed --source p --level age
#   data02_bbknn_p_age_t/*  →  data03_bbknn_p_age_t/all.txt
```

**入力例:** `data02_bbknn_r_tissue_t/Aorta.txt`, `data02_bbknn_r_tissue_t/BAT.txt`, ... (§3 の出力, `genes × samples`)
**出力例:** `data03_bbknn_r_tissue_t/all.txt` (列名は `<tissue>.txt_<sample>`)

```
@name      Skin.txt_s0    Skin.txt_s1    ...    Aorta.txt_s0    Aorta.txt_s1    ...
Sox17      0.002315       0.000415       ...    0.875226        0.827113        ...
Mrpl15     0.658357       0.660373       ...    0.433992        0.426842        ...
...
```

---

## 5. 離散化 — `05disc.py`

`04merge.py` の出力を、各遺伝子ごとに 0.1% / 75% 分位点で binary / ternary に離散化します。既定 (行方向結合) と `--transposed` (列方向結合) の双方を入力として受け付け、出力は同じディレクトリ配下に書きます。`@name` は変更しません (age 0 埋めや並び替えは行わない)。`tissue` 列の追加もしません — 入力に `tissue` 列があればそのまま素通しします (`--transposed` 入力には付きません)。定数列のみ落とします。

```sh
python 05disc.py                              # 既定 = --source r --level batch
python 05disc.py --source r --level batch     # 明示
python 05disc.py --transposed --source r --level tissue   # 列方向結合の出力を離散化
```

**入力例 (既定):** `data03_bbknn_r_batch/all.txt` (§4 既定の出力)
**入力例 (`--transposed`):** `data03_bbknn_r_tissue_t/all.txt` (§4 `--transposed` の出力, `genes × samples`)。読み込み時に `samples × genes` に転置し、もとの列ヘッダ (`<tissue>.txt_<sample>`) がそのまま `@name` 値になります。

**出力例 (既定):** `data03_bbknn_r_batch/all_disc.txt` (binary、入力の行順を保持、`tissue` 列は入力由来)

```
@name              Sox17  Mrpl15  Msc  Rdh10  ...    tissue
s0|Aorta|18m|1     1      1       1    1      ...    Aorta
s0|Aorta|24m|1     1      1       1    1      ...    Aorta
s0|Aorta|3m|1      0      1       1    1      ...    Aorta
s0|BAT|18m|1       1      1       0    1      ...    BAT
...
```

**出力例 (既定):** `data03_bbknn_r_batch/all_disc_tri.txt` (ternary, 値は 0/1/2)

```
@name              Sox17  Mrpl15  Msc  Rdh10  ...    tissue
s0|Aorta|18m|1     2      2       1    2      ...    Aorta
s0|Aorta|24m|1     2      1       1    2      ...    Aorta
...
```

**出力例 (`--transposed`):** `data03_bbknn_r_tissue_t/all_disc.txt` (`tissue` 列なし、`@name` は転置元の列名そのまま)

```
@name                       Sox17  Mrpl15  Msc  Rdh10  ...
Aorta.txt_s0|Aorta          1      0       0    0      ...
Aorta.txt_s1|Aorta          1      0       0    0      ...
BAT.txt_s0|BAT              1      1       0    1      ...
...
```

---

## 6. 離散化系バックエンド向けの整形 — `06prep_disc.py`

`05disc.py` の出力から `@name` と `tissue` を落とし、サイズ別の TSV と臓器別 TSV を切り出して `data04_bbknn_<src>_<L>{,_t}_disc/` に書き出します。FastBN / MatDNF / ilp\_bn / pgmpy の 4 つの出口がここで生成されたファイルを共通入力として使います。`--transposed` を渡すと `_t` 付きディレクトリの離散化結果を読みます (この場合 `tissue` 列がないため臓器別出力はスキップされます)。

```sh
python 06prep_disc.py                              # 既定 = --source r --level batch
python 06prep_disc.py --source r --level batch     # 明示
python 06prep_disc.py --transposed --source r --level tissue   # _t ディレクトリ
```

**入力例:** `data03_bbknn_r_batch/all_disc.txt`, `data03_bbknn_r_batch/all_disc_tri.txt` (§5 の出力)

**出力例 (すべて `data04_bbknn_r_batch_disc/` 配下):**

| ファイル                              | 内容                                |
| ------------------------------------- | ----------------------------------- |
| `all_disc.tsv`                        | 全遺伝子・binary                    |
| `all_disc{10,100,1000}.tsv`           | 先頭 N 列のサブセット               |
| `all_disc_tri.tsv` ほか               | ternary 版 (同様にサイズ別)         |
| `tissue/<tissue>.tsv`                 | 臓器別 (binary)                     |
| `tissue_tri/<tissue>.tsv`             | 臓器別 (ternary)                    |

**サンプル:** `data04_bbknn_r_batch_disc/all_disc100.tsv` (header 行 + 0/1 のみ、`@name`/`tissue` 列なし)

```
Sox17  Mrpl15  Msc  Rdh10  Stau2  Gdap1  Pi15  Il17a  Gsta3  Rims1  ...   (100 列)
1      1       1    1      1      1      1     1      1      1      ...
1      1       0    1      1      1      0     0      1      1      ...
...
```

---

## BN 推定バックエンド (5 種の出口)

### A. ingor — 大規模・連続値

`data03_bbknn_<src>_<L>_t/all.txt` (転置・連続値) を入力として、`ingor` のブートストラップ構造学習を行います。

```sh
sh 05run.sh data03_bbknn_r_tissue_t/all.txt   # §4 --transposed の出力 (resample --level tissue)
sh 05run.sh data03_bbknn_p_tissue_t/all.txt   # §4 --transposed の出力 (pseudo_bulk --level tissue)
sh 05run_each.sh                              # 各臓器ごとに個別に推定
```

**入力例:** `data03_bbknn_r_tissue_t/all.txt` (§4 `--transposed` の出力, `genes × samples`, 連続値)
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

ブートストラップ後、ECv で全体ネットワークと臓器ごとの個別ネットワークを得ます (`06run_ecv_all.sh` の `tag` 引数で `bs_all<tag>/` と `data02_bbknn<tag>_t/` を選びます。新パスを使う場合は同名のディレクトリを揃えるか `tag` を空にしておきます):

```sh
sh 06run_ecv_all.sh        # bs_all + data02_bbknn_t  → result_all.sgn3 / ecv_all/
```

**入力例:** `bs_all/result.ing.*` + `data02_bbknn_t/<tissue>.txt`
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

C++ 実装の Hill-Climb + Tabu (BIC/BDeu/K2)。`06prep_disc.py` が生成した `data04_bbknn_r_batch_disc/all_disc{,1000}.tsv` を直接読み込みます。

```sh
cd FastBN
sh compile.sh                      # g++ -O3 -std=c++17 fast_bn.cpp -o fast_bn
sh run.sh                          # → data04_bbknn_r_batch_disc/all_disc1000.tsv で構造学習
sh pred.sh                         # 学習済み構造に対する再スコア + edge-importance
```

**入力例:** `data04_bbknn_r_batch_disc/all_disc1000.tsv` (§6 の出力, 0/1 表)
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
sh run.sh                # = python learn_bn_matdnf.py ../data04_bbknn_r_batch_disc/all_disc100.tsv
```

**入力例:** `data04_bbknn_r_batch_disc/all_disc100.tsv` (§6 の出力)
**出力例:** `data04_bbknn_r_batch_disc/all_disc100_matdnf.json`

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
sh run.sh                # = python bn_ilp_pulp.py ../data04_bbknn_r_batch_disc/all_disc100.tsv
```

**入力例:** `data04_bbknn_r_batch_disc/all_disc100.tsv` (§6 の出力)
**出力例:** `data04_bbknn_r_batch_disc/all_disc100_bn_structure.json`

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
sh run.sh                # = python learn_bn_pgmpy.py ../data04_bbknn_r_batch_disc/all_disc100.tsv (hybrid+k2)
```

**入力例:** `data04_bbknn_r_batch_disc/all_disc100.tsv` (§6 の出力)
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

命名規則: `<src>` ∈ {`r` (resample), `p` (pseudo_bulk)}、`<L>` ∈ {`tissue`, `age`, `batch`}。

| ディレクトリ                        | 内容                                                    |
| ----------------------------------- | ------------------------------------------------------- |
| `data/`                             | 入力 h5ad                                               |
| `data01_<mode>/`                    | 臓器ごとの cells×genes 行列                            |
| `data02_bbknn_r_<L>/`               | resample 出力                                           |
| `data02_bbknn_p_<L>/`               | pseudo_bulk 出力                                        |
| `data02_bbknn_<src>_<L>_t/`         | 上記の転置 (genes×cells)                              |
| `data03_bbknn_<src>_<L>/all.txt`    | 行方向結合 (既定) — 離散化前の元データ                  |
| `data03_bbknn_<src>_<L>_t/all.txt`  | 列方向結合 (`--transposed`) — ingor 入口 (連続値)       |
| `data03_bbknn_<src>_<L>{,_t}/all_disc{,_tri}.txt`        | 離散化後 (binary / ternary)     |
| `data04_bbknn_<src>_<L>{,_t}_disc/all_disc{,10,100,1000}.tsv` | 離散化系 4 出口の共通入力  |
| `data04_bbknn_<src>_<L>_disc/tissue[_tri]/`              | 臓器別 0/1 (or 0/1/2) 表 (行方向結合のみ)   |
| `bs_<name>/`                        | `ingor` のブートストラップ結果                          |
| `ecv_all/`, `ecv_all2/`             | 臓器ごとの ECv 個別ネットワーク                         |
| `FastBN/`, `MatDNF/`, `ilp_bn/`, `pgmpy/` | 各バックエンドの実装                              |
