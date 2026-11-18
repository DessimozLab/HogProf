from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import polars as pl


@dataclass(frozen=True, slots=True)
class HogQueryHit:
    """One HOG similarity search hit"""

    query_hog_id: str
    hit_hog_id: str
    hit_level: str
    jaccard: float

    # Empirical p-value of `jaccard` at the hit level;
    # used for significance search only
    p_value: float | None = None


@dataclass(frozen=True, slots=True)
class HogQueryResult:
    """Search results for one or more HOG queries"""

    hits: pl.DataFrame

    @classmethod
    def from_hits(cls, hits: list[HogQueryHit], with_p_values: bool = False) -> "HogQueryResult":
        df = pl.DataFrame(
            [asdict(hit) for hit in hits],
            schema={
                "query_hog_id": pl.String,
                "hit_hog_id": pl.String,
                "hit_level": pl.String,
                "jaccard": pl.Float64,
                "p_value": pl.Float64,
            },
        )
        if not with_p_values:
            df = df.drop("p_value")
        return cls(df.sort("jaccard", descending=True))

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
    SCHEMA = {
        "hog_key": pl.Int64, "fam_id": pl.String, "hog_id": pl.String,
        "level": pl.String, "source": pl.String
    }

    def __init__(self, csv_path: Path):
        # read fam2orthoxml
        self._df = pl.read_csv(csv_path, schema_overrides=self.SCHEMA).select(self.COLUMNS)

        # check format
        if tuple(self._df.columns) != self.COLUMNS:
            raise ValueError(f"Invalid format for {csv_path}: expected {self.COLUMNS}, got {self._df.columns}")

        # Lookups are vectorized
        keys = self._df["hog_key"].to_numpy()

        # In subHOG indexes, hog_key is the row number (keys 0..n-1 in file
        # order): a key is looked up in O(1), and cannot be duplicated
        self._key_is_row = bool(np.array_equal(keys, np.arange(len(keys))))

        # Otherwise (e.g. rootHOG indexes), binary search over the sorted keys
        self._key_order = None
        self._sorted_keys = None
        if not self._key_is_row:
            self._key_order = np.argsort(keys, kind="stable")
            self._sorted_keys = keys[self._key_order]
            self._check_duplicated_keys()

    def _check_duplicated_keys(self):
        duplicated = self._sorted_keys[1:][np.diff(self._sorted_keys) == 0]
        if duplicated.size:
            raise ValueError(f"Error while loading index: duplicate hog_key {int(duplicated[0])}")

    def __len__(self) -> int:
        return self._df.height

    def _entry_at(self, i: int) -> HogEntry:
        return HogEntry(*self._df.row(int(i)))

    def _rows_where(self, column: str, value: str) -> np.ndarray:
        """Row numbers where `column` == `value`, in file order"""
        return np.flatnonzero((self._df[column] == value).to_numpy())

    def get_hog_by_key(self, key: int) -> HogEntry:
        """
        HogEntry by hog_key: O(1) if hog_key is the row number,
        binary search otherwise
        """
        if not isinstance(key, int):
            raise TypeError(f"key must be an integer, got {type(key)}")

        if self._key_is_row:
            if not 0 <= key < len(self):
                raise ValueError(f"Invalid HOG key: {key}")
            return self._entry_at(key)

        pos = int(np.searchsorted(self._sorted_keys, key))
        if pos == len(self._sorted_keys) or self._sorted_keys[pos] != key:
            raise ValueError(f"Invalid HOG key: {key}")
        return self._entry_at(self._key_order[pos])

    def get_hog_by_hog_id(self, hog_id: str) -> list[HogEntry]:
        return [self._entry_at(i) for i in self._rows_where("hog_id", hog_id)]

    def search_by_level(self, level: str) -> list[HogEntry]:
        return [self._entry_at(i) for i in self._rows_where("level", level)]

    def keys_at_level(self, level: str) -> np.ndarray:
        """hog_keys of all HOGs defined at `level`, in file order"""
        return self._df["hog_key"].to_numpy()[self._rows_where("level", level)]

    def levels(self) -> list[str]:
        """Taxonomic levels that have at least one HOG"""
        return self._df["level"].unique().sort().to_list()

    def get_level(self, level: str) -> pl.DataFrame:
        """Metadata rows of all HOGs defined at `level`, in file order"""
        return self._df.filter(pl.col("level") == level)