wget -O GSE98969_RAW.tar "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE98969&format=file"

wget -O GSE98969_experimental_design_f.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE98969&format=file&file=GSE98969%5Fexperimental%5Fdesign%5Ff%2Etxt%2Egz"

tar xvf GSE98969_RAW.tar

gunzip *.txt.gz
gunzip GSE98969_experimental_design_f.txt.gz
