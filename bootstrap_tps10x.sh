python 01preprocess_tps.py --mode 10X
python 02resample.py --target tps10x --input-glob "data01_tps10x/*.txt" ----level tissue -n 10
python 04merge.py --target tps10x --source r --level tissue
python 05disc.py --target tps10x --source r --level tissue 
python 06prep_disc.py --target tps10x --source r --level tissue


#python 03transpose.py --target tps10x --source r --level tissue
#python 04merge.py --target tps10x --source r --level tissue --transposed
#python 05disc.py --target tps10x --source r --level tissue  --transposed
#python 06prep_disc.py --target tps10x --source r --level tissue --transposed
