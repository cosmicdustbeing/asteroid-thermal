#!/usr/bin/env python3
"""
Module: observations.py

Unified reader for asteroid observation files in two formats:

1) XYZ format (per-epoch header, then lines with MJD, flux, Sun and Observer vectors):
   Header per epoch:    <num_observations> <wavelength_um>
   Lines (num_obs):     MJD  Flux  Sun_X  Sun_Y  Sun_Z  Observer_X  Observer_Y  Observer_Z
   - The header's second token is the wavelength (previously called "band").

2) RDA format (one row per observation):
   Columns:
     JD, r (AU), lambda_h (deg), beta_h (deg), delta (AU),
     lambda_o (deg), beta_o (deg), alpha (deg), wavelength_um,
     flux_mJy, flux_err_mJy
   - From these, heliocentric and observer-centric ecliptic (λ, β) are converted to
     Cartesian vectors: sun_vec (Sun→Target) and obs_vec (Observer→Target).

Outputs (normalized per observation):
  {
    'jd': float,
    'mjd': float,
    'wavelength_um': float,
    'flux_mjy': float or None,
    'flux_err_mjy': float or None,
    'sun_vec': (x, y, z),     # Sun→Target, AU
    'obs_vec': (x, y, z),     # Observer→Target, AU
    'r_au': float,            # |sun_vec|
    'delta_au': float,        # |obs_vec|
    'alpha_deg': float,       # angle between -sun_vec and obs_vec
  }

Extra utilities:
- group_sightings(gap_days=0.5): split all observations into time-contiguous groups
- tag_density(sightings, dense_min=10): label each sighting as 'dense' or 'sparse'

Notes/assumptions:
- Phase angle α is computed as angle( -sun_vec, obs_vec ).
- XYZ flux unit is treated as mJy if that’s how you store it; otherwise it is kept as given.
- RDA angles are ecliptic longitude/latitude in degrees.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import math

JD_MJD_OFFSET = 2400000.5

def distance(v: Tuple[float, float, float]) -> float:
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def dot_p(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def angle_deg(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    """Angle between vectors a and b in degrees (robust to rounding)."""
    na = distance(a)
    nb = distance(b)
    if na == 0.0 or nb == 0.0:
        return float('nan')
    c = max(-1.0, min(1.0, dot_p(a, b) / (na * nb)))
    return math.degrees(math.acos(c))

def _sph_to_cart(radius: float, lon_deg: float, lat_deg: float) -> Tuple[float, float, float]:
    """Spherical (r, λ, β) → Cartesian (x, y, z) with λ = ecliptic longitude, β = ecliptic latitude."""
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    cl = math.cos(lat)
    x = radius * cl * math.cos(lon)
    y = radius * cl * math.sin(lon)
    z = radius * math.sin(lat)
    return (x, y, z)

@dataclass
class Epoch:
    wavelength_um: float
    observations: List[Dict] = field(default_factory=list)

class Observations:
    def __init__(self, filename: str, fmt: str):
        """
        Parameters
        ----------
        filename : str
            Path to the observation file.
        fmt : str
            'xyz' or 'rda' to indicate the input format.
        """
        self.filename = filename
        self.fmt = fmt.lower().strip()
        if self.fmt not in ('xyz', 'rda'):
            raise ValueError("fmt must be 'xyz' or 'rda'")

        self.total_epochs: Optional[int] = None  # only relevant for xyz layout
        self.epochs: List[Epoch] = []           # preserves epoch grouping (xyz); rda makes 1 epoch
        self.all_observations: List[Dict] = []  # flat list across epochs for convenience

        if self.fmt == 'xyz':
            self._parse_xyz()
        else:
            self._parse_rda()

        # Compute derived quantities for all observations
        self._finalize_observations()

    # -------------------- Parsers --------------------

    def _parse_xyz(self):
        """Parse XYZ-style file with per-epoch headers and MJD-based rows."""
        try:
            with open(self.filename, 'r') as f:
                # Optional total epochs in first line
                head = f.readline()
                if not head:
                    return
                first = head.strip()
                try:
                    self.total_epochs = int(first)
                except ValueError:
                    # Not an integer: rewind and treat entire file as epoch content
                    self.total_epochs = None
                    f.seek(0)

                while True:
                    header = f.readline()
                    if not header:
                        break
                    header = header.strip()
                    if not header:
                        continue

                    parts = header.split()
                    if len(parts) < 2:
                        continue
                    try:
                        num_obs = int(parts[0])
                        wavelength_um = float(parts[1])
                    except ValueError:
                        # malformed; skip until next plausible header
                        continue

                    epoch = Epoch(wavelength_um=wavelength_um)
                    count = 0
                    while count < num_obs:
                        line = f.readline()
                        if not line:
                            break
                        line = line.strip()
                        if not line:
                            continue
                        cols = line.split()
                        if len(cols) < 8:
                            continue
                        try:
                            mjd = float(cols[0])
                            flux = float(cols[1])
                            sun_vec = (float(cols[2]), float(cols[3]), float(cols[4]))        # Sun→Target (AU)
                            obs_vec = (float(cols[5]), float(cols[6]), float(cols[7]))        # Observer→Target (AU)
                        except ValueError:
                            continue

                        obs = {
                            'jd': mjd + JD_MJD_OFFSET,
                            'mjd': mjd,
                            'wavelength_um': wavelength_um,
                            'flux_mjy': flux,
                            'flux_err_mjy': None,  # not present in xyz
                            'sun_vec': sun_vec,
                            'obs_vec': obs_vec,
                        }
                        epoch.observations.append(obs)
                        self.all_observations.append(obs)
                        count += 1

                    self.epochs.append(epoch)

        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self.filename}")

    def _parse_rda(self):
        """
        Parse RDA-style file (single table, one obs per line).

        Columns:
          0: JD
          1: r (AU)                - heliocentric distance
          2: lambda_h (deg)
          3: beta_h (deg)
          4: delta (AU)            - observer-centric distance
          5: lambda_o (deg)
          6: beta_o (deg)
          7: alpha (deg)           - provided (we will also recompute)
          8: wavelength_um
          9: flux_mJy
          10: flux_err_mJy
        """
        epoch = Epoch(wavelength_um=float('nan'))  # will vary per observation; we keep a dummy here
        try:
            with open(self.filename, 'r') as f:
                for raw in f:
                    line = raw.strip()
                    if not line or line.startswith('#'):
                        continue
                    cols = line.split()
                    if len(cols) < 11:
                        continue

                    try:
                        jd = float(cols[0])
                        r_au = float(cols[1])
                        lam_h = float(cols[2])
                        bet_h = float(cols[3])
                        delta_au = float(cols[4])
                        lam_o = float(cols[5])
                        bet_o = float(cols[6])
                        alpha_file_deg = float(cols[7])   # provided alpha (kept for reference)
                        wavelength_um = float(cols[8])
                        flux_mjy = float(cols[9])
                        flux_err_mjy = float(cols[10])
                    except ValueError:
                        continue

                    sun_vec = _sph_to_cart(r_au, lam_h, bet_h)       # Sun→Target
                    obs_vec = _sph_to_cart(delta_au, lam_o, bet_o)   # Observer→Target

                    obs = {
                        'jd': jd,
                        'mjd': jd - JD_MJD_OFFSET,
                        'wavelength_um': wavelength_um,
                        'flux_mjy': flux_mjy,
                        'flux_err_mjy': flux_err_mjy,
                        'sun_vec': sun_vec,
                        'obs_vec': obs_vec,
                        'alpha_file_deg': alpha_file_deg,  # keep the file's alpha for sanity checks
                    }
                    epoch.observations.append(obs)
                    self.all_observations.append(obs)

            # In RDA the wavelength can vary row-by-row; we leave epoch.wavelength_um as NaN
            self.epochs.append(epoch)

        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self.filename}")

    # -------------------- Post-processing --------------------

    def _finalize_observations(self):
        """Compute r, delta, and phase angle alpha for every observation."""
        for obs in self.all_observations:
            sun_vec = obs['sun_vec']
            obs_vec = obs['obs_vec']

            r_au = distance(sun_vec)
            delta_au = distance(obs_vec)
            alpha_deg = angle_deg(tuple(-c for c in sun_vec), obs_vec)

            obs['r_au'] = r_au
            obs['delta_au'] = delta_au
            obs['alpha_deg'] = alpha_deg

    # -------------------- Utilities --------------------

    def group_sightings(self, gap_days: float = 0.5) -> List[List[Dict]]:
        """
        Group observations into 'sightings' by time continuity.

        Parameters
        ----------
        gap_days : float
            Start a new group when the gap between successive observations exceeds this threshold.

        Returns
        -------
        List of lists, each inner list is a sighting (sequence of observations).
        """
        if not self.all_observations:
            return []

        obs_sorted = sorted(self.all_observations, key=lambda o: o['jd'])
        groups: List[List[Dict]] = []
        current: List[Dict] = [obs_sorted[0]]

        for prev, cur in zip(obs_sorted, obs_sorted[1:]):
            if (cur['jd'] - prev['jd']) > gap_days:
                groups.append(current)
                current = [cur]
            else:
                current.append(cur)
        groups.append(current)
        return groups

    def tag_density(self, sightings: List[List[Dict]], dense_min: int = 10) -> List[Dict]:
        """
        Tag each sighting as 'dense' or 'sparse'.

        Parameters
        ----------
        sightings : output of group_sightings()
        dense_min : int
            Minimum number of observations to be considered 'dense'.

        Returns
        -------
        List of dicts with keys: {'sighting', 'count', 'label'}
        """
        out = []
        for s in sightings:
            label = 'dense' if len(s) >= dense_min else 'sparse'
            out.append({'sighting': s, 'count': len(s), 'label': label})
        return out

    # -------------------- Representations --------------------

    def __str__(self):
        s = f"Observations from '{self.filename}' [{self.fmt.upper()}]\n"
        if self.total_epochs is not None:
            s += f"Total epochs (header): {self.total_epochs}\n"
        s += f"Epochs parsed: {len(self.epochs)}\n"
        total_obs = sum(len(ep.observations) for ep in self.epochs)
        s += f"Total observations: {total_obs}\n"
        for i, ep in enumerate(self.epochs, 1):
            w = ep.wavelength_um
            wtxt = f"{w:.3f} µm" if math.isfinite(w) else "varies"
            s += f"  Epoch {i}: {len(ep.observations)} obs, wavelength = {wtxt}\n"
        return s
