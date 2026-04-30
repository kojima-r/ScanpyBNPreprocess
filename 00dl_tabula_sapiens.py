from pathlib import Path
import re
import cellxgene_census

COLLECTION_ID = "e5f58829-1a66-40b5-a624-9046778e74f5"  # Tabula Sapiens
CENSUS_VERSION = "2025-11-08"
OUTDIR = Path("data")
OUTDIR.mkdir(exist_ok=True)

def safe_name(x):
    return re.sub(r"[^\w.-]+", "_", str(x))[:180]

with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    datasets = (
        census["census_info"]["datasets"]
        .read()
        .concat()
        .to_pandas()
    )

ts = datasets[datasets["collection_id"] == COLLECTION_ID]

print(f"Found {len(ts)} Tabula Sapiens datasets")

for _, row in ts.iterrows():
    dataset_id = row["dataset_id"]
    title = row.get("dataset_title", dataset_id)
    out = OUTDIR / f"{safe_name(title)}__{dataset_id}.h5ad"

    if out.exists():
        print("skip:", out)
        continue

    print("download:", title)
    cellxgene_census.download_source_h5ad(
        dataset_id=dataset_id,
        to_path=str(out),
        census_version=CENSUS_VERSION,
    )
