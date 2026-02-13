from single_cell_analysis import run_full_pipeline
import yaml


with open("single_cell_config.yaml", "r") as f:
    cfg = yaml.safe_load(f)

adata = run_full_pipeline(
    h5_paths=cfg["inputs"]["h5_paths"],
    sample_names=cfg["inputs"]["sample_names"],
    batch_backend=cfg["batch"]["backend"],
    out_prefix=cfg["output"]["prefix"]
)
