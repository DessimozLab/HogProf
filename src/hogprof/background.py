"""
Estimate background distribution of the Jaccard
distance per taxonomic levels. Used with --slice-subhogs
to support search for significant hits based on the
p-value threshold.
"""

from __future__ import annotations

import json
import logging
import os.path
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

from hogprof import hashstore
from hogprof.artifacts import CALIBRATION, MAPPING
from hogprof.hashstore import HashStore
from hogprof.manifest import IndexManifest
from hogprof.metadata import HogEntryIndex
from hogprof.presets import DEFAULT_ALPHAS

FORMAT_VERSION = 1

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class LevelBackground:
    # Number of indexed HOGs defined at this level
    num_hogs: int
    # Number of them sampled
    num_sampled: int
    # Null histogram: null_counts[i] pairs have similarity i / nperm
    null_counts: list[int]
    # alpha_key(alpha) -> minimal significant Jaccard, None if unreachable
    thresholds: dict[str, float | None]

    @property
    def num_pairs(self) -> int:
        return int(sum(self.null_counts))



@dataclass(frozen=True)
class Calibration:
    num_permutations: int
    sample_size: int
    seed: int
    alphas: list[float]
    levels: dict[str, LevelBackground]
    created: str = ""
    format_version: int = FORMAT_VERSION

    def threshold(self, level: str, alpha: float) -> float | None:
        """Minimal significant Jaccard at `level`, None if `alpha` cannot be reached there"""
        if alpha_key(alpha) not in {alpha_key(a) for a in self.alphas}:
            calibrated = ", ".join(str(a) for a in self.alphas)
            raise ValueError(
                f"Significance level {alpha} is not calibrated (calibrated: {calibrated}). "
                f"Run `hogprof calibrate --index-dir ... -t {alpha}`"
            )
        background = self.levels.get(level)
        if background is None:
            return None
        return background.thresholds[alpha_key(alpha)]

    def p_values(self, level: str, jaccards: np.ndarray) -> np.ndarray:
        return empirical_p_values(np.asarray(self.levels[level].null_counts),
                                  jaccards, self.num_permutations)

    def write(self, index_dir: Path) -> Path:
        path = Path(index_dir) / CALIBRATION
        data = asdict(self)
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(data) + "\n")
        tmp.replace(path)
        return path

    @classmethod
    def read(cls, path: Path) -> Calibration:
        path = Path(path)
        index_dir = path.parent
        if not path.is_file():
            raise FileNotFoundError(
                f"No {CALIBRATION} in {index_dir}. "
                f"Run `hogprof calibrate --index-dir {index_dir}` first"
            )

        data = json.loads(path.read_text())
        if data.get("format_version") != FORMAT_VERSION:
            raise ValueError(f"{path}: unsupported calibration format version "
                             f"{data.get('format_version')!r}; rerun `hogprof calibrate`")

        data["levels"] = {name: LevelBackground(**level)
                          for name, level in data["levels"].items()}
        return cls(**data)



def _validate_index(manifest: IndexManifest) -> None:
    if not manifest.slice_subhogs:
        raise ValueError("Calibration needs a subHOG index: rebuild it with --slice-subhogs")
    

def _null_histogram(packed: np.ndarray,
                    fam_ids: np.ndarray,
                    num_permutations: int) -> np.ndarray:
    """
    Makes a histogram of all-vs-all similarities for a set of packed hashes,
    excluding same-family pairs
    """
    n = len(packed)
    i, j = np.triu_indices(n, k=1)
    keep = fam_ids[i] != fam_ids[j]
    sim = hashstore.jaccard_matrix(packed)[i[keep], j[keep]]
    bins = np.rint(sim * num_permutations).astype(np.int64)
    return np.bincount(bins, minlength=num_permutations + 1)


def alpha_key(alpha: float) -> str:
    """Canonical text form of alpha, used as a JSON key and for lookups"""
    return repr(float(alpha))


def empirical_p_values(null_counts: np.ndarray, jaccards: np.ndarray,
                       num_permutations: int) -> np.ndarray:
    """
    Empirical p-values of Jaccard estimates against a null histogram
    (counts of null pairs per value 0/nperm, 1/nperm, ..., nperm/nperm)
    """
    null_counts = np.asarray(null_counts, dtype=np.int64)

    # survival[i] = number of null pairs with similarity >= i / nperm
    survival = np.cumsum(null_counts[::-1])[::-1]
    bins = np.rint(np.asarray(jaccards) * num_permutations).astype(np.int64)
    return (survival[bins] + 1) / (survival[0] + 1)


def _thresholds(null_counts: np.ndarray,
                alphas: Iterable[float],
                num_permutations: int) -> dict[str, float | None]:
    grid = np.arange(num_permutations + 1) / num_permutations
    p = empirical_p_values(null_counts, grid, num_permutations)
    result = {}
    for alpha in alphas:
        # p is non-increasing in J: the first passing value is the threshold
        passing = np.flatnonzero(p <= alpha)
        result[alpha_key(alpha)] = float(grid[passing[0]]) if passing.size else None
    return result


def calibrate(index_dir: Path,
              alphas: Iterable[float] = DEFAULT_ALPHAS,
              sample_size: int = 500,
              seed: int | None = None) -> Calibration:
    """Compute the background distributions of an index and save them next to it"""
    from hogprof.cli import track_progress

    index_dir = Path(index_dir)

    manifest = IndexManifest.read(index_dir)
    _validate_index(manifest)

    alphas = sorted({float(a) for a in alphas})
    if seed is None:
        seed = manifest.hash_seed
    rng = np.random.default_rng(seed)
    nperm = manifest.num_permutations

    hog_index = HogEntryIndex(index_dir / MAPPING)
    levels = {}
    with HashStore.from_index_dir(index_dir, manifest=manifest) as hashes:
        for level in track_progress(hog_index.levels(), description="Calibrating levels"):
            entries = hog_index.get_level(level)
            n = entries.height

            # subsample
            picked = np.sort(rng.choice(n, size=min(n, sample_size), replace=False))
            keys = entries["hog_key"].to_numpy()[picked]
            fam_ids = entries["fam_id"].to_numpy()[picked]

            # compute the null distribution
            counts = _null_histogram(hashes.get(keys), fam_ids, nperm)

            # save thresholds
            levels[level] = LevelBackground(
                num_hogs=n,
                num_sampled=len(picked),
                null_counts=counts.tolist(),
                thresholds=_thresholds(counts, alphas, nperm),
            )

    calibration = Calibration(num_permutations=nperm,
                              sample_size=sample_size,
                              seed=seed,
                              alphas=alphas,
                              levels=levels,
                              created=datetime.now(UTC).isoformat(timespec="seconds"),
    )
    path = calibration.write(index_dir)

    for alpha in alphas:
        missing = [lv for lv in levels.values() if lv.thresholds[alpha_key(alpha)] is None]
        if not missing:
            continue
        # p-values cannot go below 1 / (pairs + 1); otherwise the null is
        # too similar: even identical profiles are common among random pairs
        too_few = sum(1 / (lv.num_pairs + 1) > alpha for lv in missing)
        logger.warning("alpha=%s: no threshold at %d of %d levels "
                       "(%d with too few HOGs, %d where identical profiles are too common)",
                       alpha, len(missing), len(levels), too_few, len(missing) - too_few)
    logger.info("Calibrated %d levels at alpha %s; saved to %s",
                len(levels), ", ".join(map(str, alphas)), path)
    return calibration


def run_calibrate(index_dir: Path,
                  threshold: list[float] | None = None,
                  sample_size: int = 500,
                  seed: int | None = None) -> Calibration:
    alphas = set(DEFAULT_ALPHAS) | set(threshold or ())
    return calibrate(index_dir, alphas=alphas, sample_size=sample_size, seed=seed)
