# fitting.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
from scipy.optimize import least_squares

import shape      # for Sphere, etc.
import models     # for FRM, NEATM, ROASTM, ISO
from flux import Flux

def _pV_from_D_H(D_km: float, H_mag: float) -> float:
    # D = 1329 km * 10^(-H/5) / sqrt(pV)  -> pV = (1329*10^(-H/5)/D)^2
    if D_km <= 0:
        return np.nan
    return (1329.0 * 10.0**(-H_mag/5.0) / D_km)**2

def _phase_integral_q(G: float) -> float:
    # Bowell et al. H-G: q ≈ 0.290 + 0.684 G
    return 0.290 + 0.684*G

def _bond_albedo_from_D_HG(D_km: float, H_mag: float, G: float) -> float:
    pV = _pV_from_D_H(D_km, H_mag)
    q  = _phase_integral_q(G)
    return q * pV, pV, q

@dataclass
class FitResult:
    model: str
    shape_repr: str
    param_names: List[str]
    x_best: np.ndarray
    x_err_1sigma: Optional[np.ndarray]
    red_chisq: float
    dof: int
    success: bool
    message: str
    # optional MC
    x_mc_median: Optional[np.ndarray] = None
    x_mc_p16: Optional[np.ndarray] = None
    x_mc_p84: Optional[np.ndarray] = None

    def as_dict(self) -> Dict[str, float]:
        return {k: float(v) for k, v in zip(self.param_names, self.x_best)}

class Fitting:
    """
    Usage:
        from obs import Observations
        obs_info = Observations("2024YR4_combjwst.fluxes", fmt="rda")

        fit = Fitting(
            obs_info,
            shape_type="Sphere",
            thermal_model="FRM",
            nlat=46, nlon=92,
            model_kwargs={"Albedo": 0.1, "epsilon": 0.9}  # passed to models.<MODEL>(shape_mod, Albedo, epsilon, r_au)
        )

        result = fit.run(
            initial_params={"diam_eff_km": 1.0},
            bounds={"diam_eff_km": (0.01, 10000.0)},
            mc_trials=0
        )
    """

    def __init__(self,
                 observations: Any,
                 shape_type: str = "Sphere",
                 thermal_model: str = "FRM",
                 nlat: int = 15,
                 nlon: int = 90,
                 model_kwargs: Optional[Dict[str, Any]] = None,
                 H_mag: float | None = None,
                 G: float | None = None):
        self.obs = observations
        self.thermal_model = thermal_model.upper()
        self.shape_mod = getattr(shape, shape_type)(nlat=nlat, nlon=nlon)
        self.model_kwargs = dict(model_kwargs or {})
        # save HG (if provided we will derive Bond albedo A from D,H,G each step)
        self.H_mag = H_mag
        self.G = G

        # parameter schema per model (extend as needed)
        self.param_templates: Dict[str, List[str]] = {
            "FRM":   ["diam_eff_km"],
            "NEATM": ["diam_eff_km", "eta"],
            "ROASTM":["diam_eff_km", "eta", "f_night"],
            "ISO":   ["diam_eff_km"],
        }
        if self.thermal_model not in self.param_templates:
            raise NotImplementedError(f"Unknown model '{self.thermal_model}'")

        # normalize usable observations (must have flux and err)
        self._rows = [ob for ob in self.obs.all_observations
                      if ob.get("flux_mjy") is not None and ob.get("flux_err_mjy") is not None]
        if not self._rows:
            raise ValueError("No observations with flux and flux_err found.")

        self.F_obs = np.array([ob["flux_mjy"] for ob in self._rows], dtype=float)
        self.Sigma = np.array([max(1e-12, float(ob["flux_err_mjy"])) for ob in self._rows], dtype=float)

    # ---------- public API ----------

    def run(
        self,
        initial_params: Dict[str, float],
        bounds: Optional[Dict[str, Tuple[float, float]]] = None,
        mc_trials: int = 0,
        random_seed: Optional[int] = 12345,
    ) -> FitResult:

        pnames = self.param_templates[self.thermal_model]
        x0 = np.array([float(initial_params[p]) for p in pnames], dtype=float)

        if bounds is None:
            lo = np.full_like(x0, -np.inf)
            hi = np.full_like(x0,  np.inf)
        else:
            lo = np.array([bounds.get(p, (-np.inf, np.inf))[0] for p in pnames], dtype=float)
            hi = np.array([bounds.get(p, (-np.inf, np.inf))[1] for p in pnames], dtype=float)

        def residuals(x: np.ndarray) -> np.ndarray:
            params = {p: float(v) for p, v in zip(pnames, x)}
            F_mod = self._model_flux_vector(params)
            return (self.F_obs - F_mod) / self.Sigma

        res = least_squares(residuals, x0, bounds=(lo, hi), method="trf", jac='2-point')

        r = residuals(res.x)
        chisq = float(np.sum(r**2))
        dof = max(1, self.F_obs.size - res.x.size)
        red_chisq = chisq / dof

        # covariance estimate from J
        x_err = None
        try:
            J = res.jac
            JTJ = J.T @ J
            JTJ_inv = np.linalg.pinv(JTJ)
            s2 = chisq / dof
            cov = s2 * JTJ_inv
            x_err = np.sqrt(np.diag(cov))
        except Exception:
            x_err = None

        out = FitResult(
            model=self.thermal_model,
            shape_repr=repr(self.shape_mod),
            param_names=pnames,
            x_best=res.x.copy(),
            x_err_1sigma=(x_err.copy() if x_err is not None else None),
            red_chisq=red_chisq,
            dof=dof,
            success=bool(res.success),
            message=str(res.message),
        )

        # Monte Carlo (optional)
        if mc_trials and mc_trials > 0:
            rng = np.random.default_rng(random_seed)
            X = np.zeros((mc_trials, res.x.size), dtype=float)
            for k in range(mc_trials):
                F_pert = self.F_obs + rng.normal(0.0, self.Sigma)

                def residuals_mc(x: np.ndarray) -> np.ndarray:
                    params = {p: float(v) for p, v in zip(pnames, x)}
                    F_mod = self._model_flux_vector(params)
                    return (F_pert - F_mod) / self.Sigma

                res_mc = least_squares(residuals_mc, res.x, bounds=(lo, hi), method="trf", jac='2-point')
                X[k] = res_mc.x

            out.x_mc_median = np.median(X, axis=0)
            out.x_mc_p16 = np.percentile(X, 16, axis=0)
            out.x_mc_p84 = np.percentile(X, 84, axis=0)

        return out

    # ---------- internals ----------

    def _instantiate_model(self, r_au: float, A_bond: float):
        ModelClass = getattr(models, self.thermal_model)
        # Always pass the *current* Bond albedo derived from (D,H,G) if H,G provided;
        # otherwise fall back to user-supplied Albedo in model_kwargs (if present).
        kwargs = dict(self.model_kwargs)  # copy
        eps = self.model_kwargs.get("epsilon")
        if self.H_mag is not None and self.G is not None:
            kwargs["A"] = A_bond
        elif "A" not in kwargs:
            raise ValueError("Either provide H_mag & G (to derive Albedo) or a fixed 'Albedo' in model_kwargs.")
        if self.thermal_model == "NEATM":
            eta = params["eta"]  # could be fixed or optimized
            ss_lon = self.model_kwargs.get("ss_lon", 0.0)
            ss_lat = self.model_kwargs.get("ss_lat", 0.0)

        kwargs["r_au"] = r_au
        return ModelClass(self.shape_mod, **kwargs)

    def _model_flux_vector(self, params: dict) -> np.ndarray:
        D = float(params["diam_eff_km"])
        A_bond = None
        if self.H_mag is not None and self.G is not None:
            A_bond, pV, q = _bond_albedo_from_D_HG(D, self.H_mag, self.G)

        method_name = f"{self.thermal_model.lower()}_flux"
        F = np.empty(len(self._rows), float)

        for i, ob in enumerate(self._rows):
            lam_um = float(ob["wavelength_um"])
            Delta  = float(ob["delta_au"])
            r_au   = float(ob["r_au"])

            mdl = self._instantiate_model(r_au=r_au,
                                          A_bond=A_bond if A_bond is not None else self.model_kwargs.get("Albedo"))

            # Call the model (returns Flux object)
            flux_obj = getattr(mdl, method_name)(D, lam_um, Delta)

            if not isinstance(flux_obj, Flux):
                raise TypeError(f"{method_name} should return a Flux object, got {type(flux_obj)}")

            # Find the flux at the observation wavelength
            # wavelengths are stored in meters in Flux, convert lam_um → meters
            lam_m = lam_um * 1e-6
            # find the nearest wavelength in the flux object
            idx = np.argmin(np.abs(flux_obj.wavelengths - lam_m))
            f_w_m2_m = flux_obj.fluxes[idx]   # W/(m²·m)

            # Convert to Jy and then to mJy
            f_Jy  = flux_obj.to_Jy()[idx]     # Jy
            f_mJy = f_Jy * 1e3                # mJy

            F[i] = f_mJy

        return F
