from __future__ import annotations

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[3]
FIG_DIR = REPO_ROOT / "results_v2_corrected" / "figures_codex"
REF_DIR = REPO_ROOT / "results_v2_corrected" / "figures"


def safe_output_path(path: Path) -> Path:
    path = path.resolve()
    if "figures_codex" not in str(path) and "_codex" not in str(path):
        raise RuntimeError(f"Refusing to write outside codex isolation path: {path}")
    return path


def output_paths(stem: str) -> tuple[Path, Path]:
    pdf = safe_output_path(FIG_DIR / f"{stem}.pdf")
    png = safe_output_path(FIG_DIR / f"{stem}.png")
    return pdf, png
