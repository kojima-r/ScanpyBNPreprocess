# ScanpyBNPreprocess
Tabula Muris Senis データセットには、マウス全身のシングルセルRNA発現量データが含まれている。より具体的には、以下の２種類のデータ取得方法によるデータセットとそれらを統合し、前処理したデータセットが含まれている。
- 1：FACS（Fluorescence-Activated Cell Sorting）方式によるシングルセルRNA-seq
- 2：ドロップレット方式によるシングルセルRNA-seq
- 3：Batch Balanced k-Nearest Neighbors の略で、FACS と droplet の両データを統合した後、バッチ効果補正を行った統合データ（bbknn）



以下からデータのダウンロード:
https://figshare.com/ndownloader/files/23936555

```
mkdir data

wget --content-disposition https://ndownloader.figshare.com/files/23936555 -O data/tabula-muris-senis-bbknn-processed-official-annotations.h5ad

wget --content-disposition https://ndownloader.figshare.com/files/23937842 -O data/tabula-muris-senis-facs-processed-official-annotations.h5ad

```
#### 01. Preprocess
以下のデータを読み込み
`data/tabula-muris-senis-bbknn-processed-official-annotations.h5ad`
別データを使う場合はソースコード内を書き換えて以下を使用するようにする
`data/tabula-muris-senis-facs-processed-official-annotations.h5ad`

##### Notebook
`preprocess_facs.ipynb`
および`preprocess_bbknn.ipynb`
はh5adに標準的な処理をして、UMAPで確認したもの

### 02. resample/pseudo_bulk
(現状版ではageやbatch情報は無視している)
`resample.py`:
ブートストラップサンプリングと同様のサンプリングを行う（N=再サンプル数）

`pseudo_bulk.py`:
全細胞を平均して１サンプルにする

##### 出力

- `02data_bbknn/` :resample.pyの出力先
- `02data_bbknn2/` :pseudo_bulk.pyの出力先
  
### 03. transpose
resampleおよびpseudo_bulkの出力結果の転置を取って、
「行：遺伝子、列：サンプル」の遺伝子発現データに標準的な行列形式に変換する

##### 出力

- `02data_bbknn_t/`
- `02data_bbknn2_t/`

### 04. merge
`02data_bbknn_t/`と`02data_bbknn2_t/`
のそれぞれのファイルを読み込み、サンプル方向に結合する（全臓器のファイルを作成する、サンプル名（１行目にファイル名を付けてサンプルを区別している））

##### 出力

- `03data_bbknn/all.txt`
- `03data_bbknn/all2.txt`

### 05. BaysianNetwork推定

`03data_bbknn/all.txt`と`03data_bbknn/all2.txt`のそれぞれについて、
`05run.sh`と`05run2.sh`でBNの構造推定をおこなう

`05run_each.sh`は
`02data_bbknn_t/`を元に各臓器ごとの計算を行う
##### 出力

- `bs_all/`
- `bs_all2/`
- `bs_<file name>/`



### 06. 個別ネットワーク推定


`bs_all/result.ing.*`から全体のネットワーク`result_all.sgn3`を作成
その後、各臓器ごとの個別ネットワークを推定
`ecv_all/`に出力

##### ネットワークの可視化
- `BN_ecv.ipynb`

### XX. ADマウスとの比較
データダウンロード：
`AD/run.sh`

前処理＆個別ネットワーク計算
`AD_mouse.ipynb`


