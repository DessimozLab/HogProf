from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import polars as pl


@dataclass(frozen=True, slots=True)
class HogQueryResult:
    """Search results for one or more HOG queries"""

    hits: pl.DataFrame

    def __len__(self) -> int:
        return self.hits.height

    def to_pandas(self):
        return self.hits.to_pandas()

    def to_dicts(self) -> list[dict]:
        return self.hits.to_dicts()


@dataclass(frozen=True, slots=True)
class HogEntry:
    """
    A metadata entry for one indexed HOG. Represents a row in the fam2orthoxml file
    """
    # The integer primary key used by the LSHForest.add and hashes.h5
    hog_key: int

    # Family ID in string format:
    #   - its HOG ID for the main mode
    #   - its rootHOG's ID for --slicesubhogs
    # Example: "HOG:0015853_1"
    fam_id: str

    # HOG ID in string format:
    #   - the same as fam_id for the main mode.
    #   - subHOG ID for --slicesubhogs
    # Example: "HOG:0015853.1aa_4"
    hog_id: str

    # Taxonomic level at which the HOG is defined.
    # Example: "Archosauria"
    level: str

    # If orthoxml:
    #   - filename of the orthoxml file where the HOG was defined
    # otherwise empty
    source: str


class HogEntryIndex:
    """Indexed access to HOG metadata."""

    COLUMNS = ("hog_key", "fam_id", "hog_id", "level", "source")

    def __init__(self, csv_path: Path):
        schema = {
            "hog_key": pl.Int64, "fam_id": pl.String, "hog_id": pl.String,
            "level": pl.String, "source": pl.String
        }
        # read fam2orthoxml
        self._df = pl.read_csv(csv_path, schema_overrides=schema).select(self.COLUMNS)

        # check format
        if tuple(self._df.columns) != self.COLUMNS:
            raise ValueError(f"Invalid format for {csv_path}: expected {self.COLUMNS}, got {self._df.columns}")

        # make O(1) search indexes
        self._by_key: dict[int, int] = {}
        self._by_hog_id: dict[str, list[int]] = defaultdict(list)
        self._by_level: dict[str, list[int]] = defaultdict(list)

        for i, (key, _, hog_id, level, _) in enumerate(self._df.iter_rows()):
            if key in self._by_key:
                raise ValueError(f"Duplicate hog_key: {key}")

            self._by_key[key] = i
            self._by_hog_id[hog_id].append(i)
            self._by_level[level].append(i)

    def __len__(self) -> int:
        return self._df.height

    def _entry_at(self, i: int) -> HogEntry:
        return HogEntry(*self._df.row(i))

    def get_hog_by_key(self, key: int) -> HogEntry | None:
        i = self._by_key.get(key)
        return self._entry_at(i) if i is not None else None

    def search_by_hog_id(self, hog_id: str) -> list[HogEntry]:
        return [self._entry_at(i) for i in self._by_hog_id.get(hog_id, ())]

    def search_by_level(self, level: str) -> list[HogEntry]:
        return [self._entry_at(i) for i in self._by_level.get(level, ())]