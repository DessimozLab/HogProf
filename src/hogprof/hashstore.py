"""
Bulk access to the profile hashes stored in `hashes.h5`.

Each indexed profile is a datasketch WeightedMinHash: `num_permutations`
pairs (k, t). `build` stores a profile as one row of 2 * num_permutations
int32 values: k0, t0, k1, t1, ...

Here every (k, t) pair is packed into a single int64 by reinterpreting the
two adjacent int32 values. Two hashes agree on a permutation exactly when
their packed values are equal, so the Jaccard estimate becomes a plain
integer comparison over (n, num_permutations) arrays:

    J(a, b) = mean(a == b)

which is the same estimate as `WeightedMinHash.jaccard` (the fraction of
permutations where both k and t are equal).
"""

from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
from datasketch import WeightedMinHash

from hogprof.manifest import IndexManifest
from hogprof.utils.hashutils import HASH_SEED

# dtype of the values in hashes.h5
STORED_DTYPE = np.int32


def pack(rows: np.ndarray) -> np.ndarray:
    """
    Pack stored rows of shape (..., 2 * nperm) into (..., nperm) int64,
    one value per (k, t) pair
    """
    rows = np.ascontiguousarray(rows, dtype=STORED_DTYPE)
    if rows.shape[-1] % 2:
        raise ValueError(f"Row width must be even, got {rows.shape[-1]}")
    return rows.view(np.int64)


def unpack(packed: np.ndarray) -> np.ndarray:
    """Inverse of `pack`: (..., nperm) int64 -> (..., nperm, 2) int64 (k, t) pairs"""
    packed = np.ascontiguousarray(packed, dtype=np.int64)
    pairs = packed.view(STORED_DTYPE).reshape(*packed.shape, 2)
    return pairs.astype(np.int64)


def to_minhash(packed: np.ndarray, seed: int = HASH_SEED) -> WeightedMinHash:
    """Make a datasketch WeightedMinHash from one packed hash (nperm,)"""
    packed = np.asarray(packed)
    if packed.ndim != 1:
        raise ValueError(f"Expected one packed hash of shape (nperm,), got {packed.shape}")
    return WeightedMinHash(seed=seed, hashvalues=unpack(packed))


def from_minhash(minhash: WeightedMinHash) -> np.ndarray:
    """Pack a datasketch WeightedMinHash the same way `build` stores it"""
    values = np.asarray(minhash.hashvalues)
    stored = values.astype(STORED_DTYPE)
    if not np.array_equal(stored, values):
        raise ValueError("Hash values do not fit the stored int32 format")
    return pack(stored.ravel())


def jaccard(query: np.ndarray, hashes: np.ndarray) -> np.ndarray:
    """
    Jaccard estimates of one packed hash (nperm,) against
    packed hashes (n, nperm). Returns (n,) float64
    """
    query = np.asarray(query)
    hashes = np.asarray(hashes)
    if hashes.ndim == 1:
        hashes = hashes[None, :]
    if query.ndim != 1 or hashes.shape[1] != query.shape[0]:
        raise ValueError(f"Shape mismatch: query {query.shape}, hashes {hashes.shape}")
    return np.count_nonzero(hashes == query, axis=1) / query.shape[0]


def jaccard_matrix(a: np.ndarray, b: np.ndarray | None = None,
                   *, chunk_bytes: int = 64 * 2**20) -> np.ndarray:
    """
    All-vs-all Jaccard estimates between packed hashes a (n, nperm) and
    b (m, nperm), or a against itself if b is None. Returns (n, m) float64.

    Computed in row blocks, so memory stays at about `chunk_bytes`
    on top of the result.
    """
    a = np.asarray(a)
    b = a if b is None else np.asarray(b)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[1]:
        raise ValueError(f"Shape mismatch: {a.shape} vs {b.shape}")

    n, nperm = a.shape
    m = b.shape[0]
    out = np.empty((n, m), dtype=np.float64)
    # a (rows, 1, nperm) == b (1, m, nperm) takes rows * m * nperm bytes
    rows = max(1, chunk_bytes // max(1, m * nperm))
    for start in range(0, n, rows):
        block = a[start:start + rows, None, :] == b[None, :, :]
        out[start:start + rows] = np.count_nonzero(block, axis=2) / nperm
    return out


class HashStore:
    """
    Read-only access to the hashes of an index, addressed by `hog_key`
    (the same key as in fam2orthoxml.csv and the LSH forest).

    Note that hashes.h5 has unused zero rows: it grows in chunks, and in the
    rootHOG mode keys of skipped families are never written. Use keys from
    the HOG metadata, not row numbers.
    """

    def __init__(self, path: Path, *,
                 num_permutations: int,
                 dataset: str | None = None,
                 seed: int = HASH_SEED,
                 preload: bool = False):
        """
        path: hashes.h5
        num_permutations: number of permutations the index was built with
        dataset: dataset name in the file; the only dataset if None
        seed: seed of the generator the index was built with
        preload: read all hashes into memory. Faster when many hashes
            are read (scans, calibration), costs 8 * nperm bytes per row
        """
        self.path = Path(path)
        self.num_permutations = num_permutations
        self.seed = seed

        self._file = h5py.File(self.path, mode="r")
        try:
            if dataset is None:
                names = list(self._file.keys())
                if len(names) != 1:
                    raise ValueError(f"{self.path}: expected one dataset, found {names}; "
                                     f"specify which one to use")
                dataset = names[0]
            if dataset not in self._file:
                raise ValueError(f"{self.path}: no dataset {dataset!r}")
            self._dataset = self._file[dataset]

            width = self._dataset.shape[1]
            if width != 2 * num_permutations:
                raise ValueError(
                    f"{self.path}: rows have {width} values, expected "
                    f"{2 * num_permutations} for {num_permutations} permutations"
                )
            self.dataset_name = dataset
            self._packed = pack(self._dataset[:]) if preload else None
        except BaseException:
            self._file.close()
            raise

    @classmethod
    def from_index_dir(cls, index_dir: Path, *, manifest: IndexManifest | None = None,
                       preload: bool = False) -> HashStore:
        """Open the hashes of an index built by `hogprof build`"""
        index_dir = Path(index_dir)
        if manifest is None:
            manifest = IndexManifest.read(index_dir)
        return cls(index_dir / "hashes.h5",
                   num_permutations=manifest.num_permutations,
                   dataset=manifest.hashes_dataset,
                   seed=manifest.hash_seed,
                   preload=preload)

    @property
    def num_rows(self) -> int:
        """Number of rows in the dataset, including unused ones"""
        return self._dataset.shape[0]

    @property
    def preloaded(self) -> bool:
        return self._packed is not None

    def _check_keys(self, keys: np.ndarray) -> None:
        if keys.size and (keys.min() < 0 or keys.max() >= self.num_rows):
            bad = keys[(keys < 0) | (keys >= self.num_rows)]
            raise KeyError(f"HOG keys out of range [0, {self.num_rows}): {bad[:5].tolist()}")

    def get(self, keys) -> np.ndarray:
        """
        Packed hashes for `keys` (any order, duplicates allowed),
        shape (len(keys), nperm) int64
        """
        keys = np.asarray(keys)
        if keys.ndim != 1:
            raise ValueError(f"keys must be one-dimensional, got shape {keys.shape}")
        if keys.size and not np.issubdtype(keys.dtype, np.integer):
            raise TypeError(f"keys must be integers, got {keys.dtype}")
        keys = keys.astype(np.int64, copy=False)
        self._check_keys(keys)

        if self._packed is not None:
            return self._packed[keys]

        if keys.size == 0:
            return np.empty((0, self.num_permutations), dtype=np.int64)

        # h5py fancy indexing needs unique increasing indices
        unique, inverse = np.unique(keys, return_inverse=True)
        return pack(self._dataset[unique, :])[inverse]

    def get_one(self, key: int) -> np.ndarray:
        """Packed hash of one HOG, shape (nperm,)"""
        return self.get([int(key)])[0]

    def get_minhash(self, key: int) -> WeightedMinHash:
        """The hash of one HOG as a datasketch WeightedMinHash (e.g. for LSH queries)"""
        return to_minhash(self.get_one(key), seed=self.seed)

    def jaccard(self, query_key: int, keys) -> np.ndarray:
        """Jaccard estimates between one HOG and many, shape (len(keys),)"""
        return jaccard(self.get_one(query_key), self.get(keys))

    def close(self) -> None:
        self._file.close()

    def __enter__(self) -> HashStore:
        return self

    def __exit__(self, *exc) -> None:
        self.close()
