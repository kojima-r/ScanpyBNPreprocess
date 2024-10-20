# ScanpyBNPreprocess

#### 01. Preprocess
以下のデータを読み込み
`tabula-muris-senis-bbknn-processed-official-annotations.h5ad`
別データを使う場合はソースコード内を書き換えて以下を使用するようにする
`tabula-muris-senis-facs-processed-official-annotations.h5ad`


### 02. resample/pseudo_bulk

出力

- `02data_bbknn/` resample
- `02data_bbknn2/` pseudo_bulk
  
### 03. transpose

出力

- `02data_bbknn_t/`
- `02data_bbknn2_t/`

### 04. merge

出力

- `03data_bbknn/all.txt`
- `03data_bbknn/all2.txt`

### 05. BN

data=all
path=03data_bbknn
