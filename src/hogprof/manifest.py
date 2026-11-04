"""
The index manifest: a JSON file written by `build` next to the other index
files. It records the parameters needed to read the index back (so `search`
doesn't have to be told them again) and how the index was built.

The manifest is written last, so its presence also marks a completed build.
"""

from __future__ import annotations

import json
from dataclasses import MISSING, asdict, dataclass, field, fields
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

MANIFEST_NAME = "manifest.json"

# Bump when the index layout changes in a way readers must know about
FORMAT_VERSION = 1


@dataclass(frozen=True)
class IndexManifest:
    # Number of hash functions per profile; each hash is a (k, t) pair,
    # so a row in hashes.h5 has 2 * num_permutations values
    num_permutations: int

    # Seed of the WeightedMinHashGenerator used to make the hashes
    hash_seed: int

    # Name of the dataset in hashes.h5 that holds the hashes
    hashes_dataset: str

    # True if the index has a profile per subHOG (--slice-subhogs),
    # False if it has a profile per rootHOG
    slice_subhogs: bool

    # Number of indexed profiles (rows in fam2orthoxml.csv)
    num_hogs: int

    # Build options and inputs, for provenance only
    build: dict[str, Any] = field(default_factory=dict)

    hogprof_version: str = ""
    created: str = ""
    format_version: int = FORMAT_VERSION

    def write(self, index_dir: Path) -> Path:
        path = Path(index_dir) / MANIFEST_NAME
        data = asdict(self)
        if not data["created"]:
            data["created"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        if not data["hogprof_version"]:
            from hogprof import __version__
            data["hogprof_version"] = __version__

        # write-then-rename, so a crash never leaves a half-written manifest
        tmp = path.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(data, indent=2, default=str) + "\n")
        tmp.replace(path)
        return path

    @classmethod
    def read(cls, index_dir: Path) -> IndexManifest:
        path = Path(index_dir) / MANIFEST_NAME
        if not path.is_file():
            raise FileNotFoundError(
                f"No {MANIFEST_NAME} in {index_dir}. The index is incomplete or "
                f"was built by an older HogProf version; rebuild it with `hogprof build`"
            )

        data = json.loads(path.read_text())
        version = data.get("format_version")
        if version != FORMAT_VERSION:
            raise ValueError(
                f"{path}: unsupported index format version {version!r} "
                f"(this HogProf reads version {FORMAT_VERSION}); rebuild the index"
            )

        missing = [f.name for f in fields(cls)
                   if f.name not in data
                   and f.default is MISSING and f.default_factory is MISSING]
        if missing:
            raise ValueError(f"{path}: missing required fields: {', '.join(missing)}")

        # ignore unknown keys, so that newer optional fields don't break older readers
        return cls(**{f.name: data[f.name] for f in fields(cls) if f.name in data})

    @staticmethod
    def remove(index_dir: Path) -> None:
        """Remove a stale manifest before rebuilding into the same directory"""
        (Path(index_dir) / MANIFEST_NAME).unlink(missing_ok=True)

