"""Shared paths, constants, grid definition and I/O helpers for the Lagos heat project."""

from __future__ import annotations

import io
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw"
INTERIM = ROOT / "data" / "interim"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"
WORK = RAW / "_work"  # scratch arrays (gitignored with the rest of data/raw)
SUMMARIES = INTERIM / "summaries"

PC_STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
LAGOS_BBOX = (2.70, 6.35, 4.40, 6.75)  # lon/lat, generous box around Lagos State
UTM31N = "EPSG:32631"
CELL = 500  # grid cell size, metres

# Analysis periods: pooled dry seasons (December to February) per period.
PERIODS = {
    "early": ["2016-12-01/2017-02-28", "2017-12-01/2018-02-28"],
    "late": ["2023-12-01/2024-02-29", "2024-12-01/2025-02-28"],
}
MIN_OBS = 3  # cells with fewer valid observations per pixel (on average) are flagged


# ------------------------------------------------------------------ grid
@dataclass(frozen=True)
class GridSpec:
    """Regular 500 m grid in UTM 31N. Cell id = row * ncol + col, row 0 at the top."""

    xmin: float
    ymax: float
    ncol: int
    nrow: int
    cell: float = CELL

    @property
    def ncells(self):
        return self.ncol * self.nrow

    @property
    def bounds(self):
        return (self.xmin, self.ymax - self.nrow * self.cell,
                self.xmin + self.ncol * self.cell, self.ymax)

    def save(self, path):
        Path(path).write_text(json.dumps(self.__dict__, indent=2))

    @classmethod
    def load(cls, path=None):
        path = path or INTERIM / "grid_spec.json"
        return cls(**json.loads(Path(path).read_text()))

    def cell_of_xy(self, x, y):
        """Cell id for UTM coordinates (arrays); -1 outside the grid."""
        col = np.floor((np.asarray(x) - self.xmin) / self.cell).astype(np.int64)
        row = np.floor((self.ymax - np.asarray(y)) / self.cell).astype(np.int64)
        inside = (col >= 0) & (col < self.ncol) & (row >= 0) & (row < self.nrow)
        return np.where(inside, row * self.ncol + col, -1)

    def pixel_cells(self, transform, shape):
        """Cell id of each pixel centre for a north-up raster in UTM 31N."""
        nrows, ncols = shape
        xs = transform.c + (np.arange(ncols) + 0.5) * transform.a
        ys = transform.f + (np.arange(nrows) + 0.5) * transform.e
        col = np.floor((xs - self.xmin) / self.cell).astype(np.int64)
        row = np.floor((self.ymax - ys) / self.cell).astype(np.int64)
        col = np.where((col >= 0) & (col < self.ncol), col, -1).astype(np.int32)
        row = np.where((row >= 0) & (row < self.nrow), row, -1).astype(np.int32)
        cells = row[:, None] * np.int32(self.ncol) + col[None, :]
        cells[(row[:, None] < 0) | (col[None, :] < 0)] = -1
        return cells

    def pixel_cells_reprojected(self, transform, shape, crs):
        """Cell id of each pixel centre for a raster in another CRS (pixel centres reprojected)."""
        from pyproj import Transformer

        nrows, ncols = shape
        jj, ii = np.meshgrid(np.arange(ncols) + 0.5, np.arange(nrows) + 0.5)
        xs, ys = transform * (jj, ii)
        tr = Transformer.from_crs(crs, UTM31N, always_xy=True)
        ux, uy = tr.transform(xs, ys)
        return self.cell_of_xy(ux, uy)

    def geobox(self, resolution, anchor=0.0):
        """odc GeoBox covering the grid; pixel edges at multiples of `resolution` + anchor."""
        from affine import Affine
        from odc.geo.geobox import GeoBox

        x0 = np.floor((self.xmin - anchor) / resolution) * resolution + anchor
        y1 = np.ceil((self.ymax - anchor) / resolution) * resolution + anchor
        xmax, ymin = self.bounds[2], self.bounds[1]
        nx = int(np.ceil((xmax - x0) / resolution))
        ny = int(np.ceil((y1 - ymin) / resolution))
        return GeoBox((ny, nx), Affine(resolution, 0, x0, 0, -resolution, y1), UTM31N)


def aggregate(cells, values=None, mask=None, ncells=None):
    """Per-cell (sum, count) of `values` over pixels with cells >= 0 and mask True.

    `cells` gives the cell id of each pixel (-1 = outside). With values=None the sum
    is the count of masked pixels.
    """
    cells = np.asarray(cells).ravel()
    keep = cells >= 0
    if mask is not None:
        keep &= np.asarray(mask).ravel()
    idx = cells[keep]
    count = np.bincount(idx, minlength=ncells).astype("float64")
    if values is None:
        return count, count
    total = np.bincount(idx, weights=np.asarray(values, dtype="float64").ravel()[keep],
                        minlength=ncells)
    return total, count


def safe_div(a, b):
    a, b = np.asarray(a, dtype="float64"), np.asarray(b, dtype="float64")
    out = np.full(np.broadcast(a, b).shape, np.nan)
    np.divide(a, b, out=out, where=b > 0)
    return out


# --------------------------------------------------------------- outputs
def write_summary(name, data):
    SUMMARIES.mkdir(parents=True, exist_ok=True)
    path = SUMMARIES / f"{name}.json"
    path.write_text(json.dumps(data, indent=2, default=_json_default), encoding="utf-8")
    return path


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return o.as_posix()
    return str(o)


def load_grid():
    import geopandas as gpd

    return gpd.read_file(INTERIM / "grid.gpkg")


def write_raster(path, array, transform, crs=UTM31N, nodata=np.nan):
    import rasterio

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff", "height": array.shape[0], "width": array.shape[1], "count": 1,
        "dtype": array.dtype.name, "crs": crs, "transform": transform, "nodata": nodata,
        "compress": "deflate", "tiled": True, "blockxsize": 512, "blockysize": 512,
    }
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(array, 1)


# ------------------------------------------------------------ day stacks
class DayStack:
    """Disk-backed stack of per-day rasters (int16, scaled) with strip-wise medians.

    Keeps memory low on large grids: each day is written once, and the median is
    computed over row strips. Progress is recorded in a JSON sidecar so an
    interrupted run resumes where it stopped.
    """

    NODATA = -32768

    def __init__(self, name, keys, shape, scale):
        WORK.mkdir(parents=True, exist_ok=True)
        self.path = WORK / f"{name}.int16"
        self.meta_path = WORK / f"{name}.json"
        self.keys, self.shape, self.scale = list(keys), tuple(shape), scale
        meta = json.loads(self.meta_path.read_text()) if self.meta_path.exists() else None
        if (meta and meta["keys"] != self.keys and tuple(meta["shape"]) == self.shape
                and set(meta["keys"]) <= set(self.keys) and self.path.exists()):
            meta = self._extend(meta)
        fresh = not (meta and meta["keys"] == self.keys and tuple(meta["shape"]) == self.shape
                     and self.path.exists())
        self.done = {} if fresh else meta["done"]
        self.arr = np.memmap(self.path, dtype="int16", mode="w+" if fresh else "r+",
                             shape=(max(len(self.keys), 1), *self.shape))
        self._save_meta()

    def _extend(self, meta):
        """Copy stored days into a stack with additional keys (no stored day is lost)."""
        old = np.memmap(self.path, dtype="int16", mode="r",
                        shape=(max(len(meta["keys"]), 1), *self.shape))
        tmp = self.path.with_suffix(".extending")
        new = np.memmap(tmp, dtype="int16", mode="w+", shape=(len(self.keys), *self.shape))
        for k in meta["done"]:
            new[self.keys.index(k)] = old[meta["keys"].index(k)]
        new.flush()
        del old, new
        tmp.replace(self.path)
        meta = {"keys": self.keys, "shape": list(self.shape), "done": meta["done"]}
        self.meta_path.write_text(json.dumps(meta))
        return meta

    def _save_meta(self):
        self.meta_path.write_text(json.dumps(
            {"keys": self.keys, "shape": self.shape, "done": self.done}))

    def is_done(self, key):
        return key in self.done

    def forget(self, key):
        self.done.pop(key, None)
        self._save_meta()

    def write(self, key, values, info=None):
        """Store a float array (NaN = no data) for `key`."""
        out = np.full(self.shape, self.NODATA, dtype="int16")
        ok = np.isfinite(values)
        out[ok] = np.clip(np.round(values[ok] * self.scale), -32767, 32767)
        self.arr[self.keys.index(key)] = out
        self.arr.flush()
        self.done[key] = info or {}
        self._save_meta()

    def _valid(self, vals, vmin, vmax):
        ok = vals != self.NODATA
        if vmin is not None:
            ok &= vals >= vmin * self.scale
        if vmax is not None:
            ok &= vals <= vmax * self.scale
        return ok

    def slice_medians(self, keys, mask, vmin=None, vmax=None):
        """{key: median of stored values over pixels where mask is True} (None if empty).

        Values outside [vmin, vmax] are treated as missing.
        """
        out = {}
        for k in keys:
            vals = self.arr[self.keys.index(k)][mask]
            vals = vals[self._valid(vals, vmin, vmax)]
            out[k] = float(np.median(vals)) / self.scale if vals.size else None
        return out

    def count_outside(self, keys, mask, vmin, vmax):
        """{key: (n below vmin, n above vmax)} over valid stored values where mask is True."""
        out = {}
        for k in keys:
            vals = self.arr[self.keys.index(k)][mask]
            vals = vals[vals != self.NODATA]
            out[k] = (int((vals < vmin * self.scale).sum()), int((vals > vmax * self.scale).sum()))
        return out

    def median(self, keys=None, strip=64, offsets=None, vmin=None, vmax=None, q=None):
        """Per-pixel median (or percentile `q`, 0-100) over `keys` (default all) and count of
        valid days.

        Values outside [vmin, vmax] are treated as missing. With `offsets` ({key: value}),
        each day's value minus its offset is used (after the range test).
        """
        import warnings

        keys = list(keys or self.keys)
        idx = [self.keys.index(k) for k in keys]
        med = np.full(self.shape, np.nan, dtype="float32")
        cnt = np.zeros(self.shape, dtype="uint16")
        if not idx:
            return med, cnt
        off = None
        if offsets is not None:
            off = np.array([np.nan if offsets[k] is None else offsets[k] * self.scale
                            for k in keys], dtype="float32")[:, None, None]
        for r0 in range(0, self.shape[0], strip):
            raw = self.arr[idx, r0:r0 + strip, :]
            block = raw.astype("float32")
            block[~self._valid(raw, vmin, vmax)] = np.nan
            del raw
            if off is not None:
                block -= off
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                stat = (np.nanmedian(block, axis=0) if q is None
                        else np.nanpercentile(block, q, axis=0))
                med[r0:r0 + strip] = stat / self.scale
            cnt[r0:r0 + strip] = np.isfinite(block).sum(axis=0)
        return med, cnt


# ------------------------------------------------------------------ http
def retry(fn, *args, attempts=8, base=10, cap=300, **kwargs):
    """Call fn, retrying on any exception with exponential backoff (about 20 min in total).

    Used around network reads so that short outages do not end long runs or, worse,
    leave silently incomplete results.
    """
    import time

    for k in range(attempts):
        try:
            return fn(*args, **kwargs)
        except Exception as exc:  # noqa: BLE001
            if k == attempts - 1:
                raise
            wait = min(base * 2**k, cap)
            print(f"    retry {k + 1}/{attempts - 1} in {wait} s after {type(exc).__name__}: "
                  f"{str(exc)[:120]}", flush=True)
            time.sleep(wait)


def fetch(url, dest=None, params=None, timeout=300):
    """GET a URL with requests (retried); return bytes, or stream to dest and return its path."""
    return retry(_fetch, url, dest, params, timeout)


def _fetch(url, dest=None, params=None, timeout=300):
    import requests

    with requests.get(url, params=params, timeout=timeout, stream=dest is not None) as r:
        r.raise_for_status()
        if dest is None:
            return r.content
        dest = Path(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        tmp = dest.with_suffix(dest.suffix + ".part")
        with open(tmp, "wb") as f:
            for chunk in r.iter_content(1 << 20):
                f.write(chunk)
        tmp.replace(dest)
        return dest


class HTTPRangeFile(io.RawIOBase):
    """Read-only, seekable view of a remote file using HTTP Range requests.

    Lets rasterio read windows of a remote GeoTIFF through requests instead of
    GDAL's curl, whose Windows (schannel) build fails certificate revocation
    checks for some hosts. Reads are served from 512 KiB blocks cached per URL.
    Only works if the server honours Range.
    """

    BLOCK = 1 << 19
    session = None
    _sizes: dict = {}
    _cache: dict = {}
    bytes_fetched = 0

    def __init__(self, url):
        import requests

        self.url, self.pos = url, 0
        if HTTPRangeFile.session is None:
            HTTPRangeFile.session = requests.Session()
        if url not in self._sizes:
            r = self.session.head(url, allow_redirects=True, timeout=60)
            r.raise_for_status()
            self._sizes[url] = int(r.headers["Content-Length"])
        self.size = self._sizes[url]

    @classmethod
    def clear_cache(cls):
        cls._cache.clear()

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def _block(self, i):
        key = (self.url, i)
        if key not in self._cache:
            start = i * self.BLOCK
            end = min(start + self.BLOCK, self.size) - 1
            for attempt in range(4):
                try:
                    r = self.session.get(self.url, headers={"Range": f"bytes={start}-{end}"},
                                         timeout=120)
                    break
                except Exception:  # noqa: BLE001
                    if attempt == 3:
                        raise
            if r.status_code != 206:
                raise OSError(f"Range request not honoured (HTTP {r.status_code}): {self.url}")
            self._cache[key] = r.content
            HTTPRangeFile.bytes_fetched += len(r.content)
        return self._cache[key]

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        n = min(n, self.size - self.pos)
        if n <= 0:
            return b""
        out, p, end = bytearray(), self.pos, self.pos + n
        while p < end:
            i = p // self.BLOCK
            chunk = self._block(i)[p - i * self.BLOCK : p - i * self.BLOCK + (end - p)]
            out += chunk
            p += len(chunk)
        self.pos = p
        return bytes(out)

    def readinto(self, b):
        data = self.read(len(b))
        b[: len(data)] = data
        return len(data)


def http_opener(path, mode="rb"):
    """rasterio `opener` for https URLs. rasterio probes with a dummy path first."""
    if not str(path).startswith("http"):
        raise FileNotFoundError(path)
    return HTTPRangeFile(path)


def pc_client():
    import planetary_computer
    import pystac_client

    return pystac_client.Client.open(PC_STAC, modifier=planetary_computer.sign_inplace)


def resign(item):
    """Copy of a Planetary Computer item with fresh SAS tokens. Tokens expire after about an
    hour, so long runs re-sign before every load attempt instead of reusing search-time URLs."""
    import planetary_computer

    item = item.clone()
    for asset in item.assets.values():
        asset.href = asset.href.split("?", 1)[0]
    return planetary_computer.sign(item)


def configure_gdal():
    """GDAL settings for reading cloud-optimised GeoTIFFs over HTTP."""
    import odc.stac

    odc.stac.configure_rio(cloud_defaults=True, GDAL_HTTP_MAX_RETRY="5",
                           GDAL_HTTP_RETRY_DELAY="2")


def lagos_lgas():
    """Full-resolution Lagos LGA polygons (written by grid.py), in UTM 31N."""
    import geopandas as gpd

    return gpd.read_file(INTERIM / "lagos_lgas.gpkg")
