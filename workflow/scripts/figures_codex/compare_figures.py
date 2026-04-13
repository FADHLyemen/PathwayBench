from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PIL import Image
from skimage.metrics import structural_similarity as ssim

from common_figures import FIG_DIR, REF_DIR, safe_output_path


def compute_ssim(fig_num: int) -> float:
    a = Image.open(FIG_DIR / f"Figure_{fig_num}.png").convert("L")
    b = Image.open(REF_DIR / f"Figure_{fig_num}.png").convert("L")
    if a.size != b.size:
      b = b.resize(a.size)
    arr_a = np.asarray(a)
    arr_b = np.asarray(b)
    return float(ssim(arr_a, arr_b, data_range=255))


def main():
    out = {f"Figure_{n}": compute_ssim(n) for n in [1, 2, 3, 4, 5, 6, 8]}
    target = safe_output_path(FIG_DIR / "ssim_results.json")
    target.write_text(json.dumps(out, indent=2))
    for k, v in out.items():
        print(f"{k}\t{v:.4f}")


if __name__ == "__main__":
    main()
