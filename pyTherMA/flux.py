#!/usr/bin/env python

# flux.py
"""
This module provides a function to compute the disk-integrated flux
from a shape model that has facet temperatures, areas, and emission angles.
This flux computation is decoupled from the thermal model so that
the same routine can be used regardless of which thermal model is applied.
"""

import numpy as np
from constants import AU, H, C, K_B

def planck_radiance(wavelength, T):
    """
    Compute the spectral radiance (B_lambda) using the Planck function.
    Parameters:
        wavelength (float or numpy.ndarray): Wavelength in micrometers.
        T (float or numpy.ndarray): Temperature in Kelvin.
    Returns:
        numpy.ndarray: Spectral radiance in W/(m²·m·sr)
    """
    # Ensure wavelength and T are numpy arrays (and convert wavelength to meters).
    wavelength = np.asarray(wavelength)*1e-6
    T = np.asarray(T)

    T_adj = np.where(T == 0, 2.7, T)

    # Calculate the exponent
    #exponent = (H * C) / (wavelength * K_B * T)
    exponent = np.where(T_adj > 0, (H * C) / (wavelength * K_B * T_adj), np.inf)
    #exponent = np.divide(H * C, wavelength * K_B * T,
    #                 out=np.full_like(T, np.inf), where=T > 0)
    # Use numpy exp safely (for T==0, set radiance to 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        B_lambda = (2 * H * C**2) / (wavelength**5) / (np.exp(exponent) - 1)
    B_lambda = np.where(T > 0, B_lambda, 0)
    return B_lambda

def flux_to_Jy(flux, wavelength):
    """
    Convert spectral flux density from W/(m²·m) to Jansky.
    
    Parameters:
        flux (float or np.ndarray): Spectral flux density in W/(m²·m).
        wavelength (float or np.ndarray): Wavelength in meters at which flux is measured.
        
    Returns:
        np.ndarray: Spectral flux density in Jansky (Jy), where 1 Jy = 1e-26 W/(m²·Hz).
    
    The conversion is done using the relation:
        F(ν) = F(λ) * (λ² / c)
        F_Jy = F(ν) / 1e-26
    """
    flux = np.asarray(flux)

    # Convert wavelength to meter.
    wavelength_m = np.asarray(wavelength)*1e-6
    
    # Convert flux density from per meter to per Hz
    F_nu = flux * (wavelength_m**2 / C)
    
    # Convert to Jansky:
    flux_Jy = F_nu / 1e-26
    return flux_Jy

class Flux:
    """
    A class to encapsulate flux computations.
    """
    def __init__(self, wavelengths, fluxes):
        """
        Initialize the Flux object.
        
        Parameters:
            wavelengths (array-like): Wavelengths in meters.
            fluxes (array-like): Flux densities in W/(m²·m).
        """
        self.wavelengths = np.asarray(wavelengths)
        self.fluxes = np.asarray(fluxes)
    
    def to_Jy(self):
        """
        Convert the stored fluxes to Jansky.
        
        Returns:
            np.ndarray: Fluxes in Jansky.
        """
        return flux_to_Jy(self.fluxes, self.wavelengths)
    
    def __str__(self):
        s = "Flux Object:\n"
        s += f"  Wavelengths (micrometers): {self.wavelengths}\n"
        s += f"  Fluxes (W/(m²·m)): {self.fluxes}\n"
        return s

class Lightcurve:
    """
    A class to store lightcurve fluxes.
    """
    def __init__(self, rotation_phases, wavelengths, fluxes):
        """
        Initialize the Lightcurve object.

        Parameters:
            rotation_phases (array_like): Rotation phases in degrees.
            wavelengths (array-like): Wavelengths in meters.
            fluxes (array-like): Flux densities in W/(m²·m).
        """
        rp = np.asarray(rotation_phases, dtype=float).ravel()
        wl = np.asarray(wavelengths, dtype=float).ravel()
        Fx = np.asarray(fluxes, dtype=float)

        # Normalize flux shape to (N, M)
        if Fx.ndim == 1:
            # allow (N,) only if exactly one wavelength
            if wl.size != 1 or Fx.size != rp.size:
                raise ValueError("1D fluxes require exactly one wavelength and len(fluxes)==len(rotation_phases).")
            Fx = Fx.reshape(rp.size, 1)
        elif Fx.ndim == 2:
            N, M = Fx.shape
            if N != rp.size:
                raise ValueError(f"fluxes.shape[0] ({N}) must match number of rotation phases ({rp.size}).")
            if M != wl.size:
                raise ValueError(f"fluxes.shape[1] ({M}) must match number of wavelengths ({wl.size}).")
        else:
            raise ValueError("fluxes must be 1D or 2D (N) or (N, M).")

        self.rotation_phases = rp
        self.wavelengths = wl
        self.fluxes = Fx  # (N, M)

    def to_Jy(self):
        """Return fluxes converted to Jansky with same (N, M) shape."""
        n_phase, n_lambda = self.fluxes.shape
        fluxes_Jy = np.zeros_like(self.fluxes, dtype=float)

        for j in range(n_lambda):
            # flux_to_Jy expects flux (array of shape (N,)) and wavelength (scalar or array in µm)
            fluxes_Jy[:, j] = flux_to_Jy(self.fluxes[:, j], self.wavelengths[j])

        return fluxes_Jy
            
    def select_wavelength(self, wavelength, atol=0.0):
        """
        Return (phases, flux_at_lambda) for the closest stored wavelength.
        If atol>0, require |Δλ| <= atol.
        """
        idx = int(np.argmin(np.abs(self.wavelengths - wavelength)))
        if atol > 0 and np.abs(self.wavelengths[idx] - wavelength) > atol:
            raise ValueError("Requested wavelength not found.")
        return self.rotation_phases.copy(), self.fluxes[:, idx].copy()

    def at_phase(self, phase_deg):
        """
        Get a Flux object slice at a given rotation phase (nearest neighbor).
        Returns a Flux(wavelengths, fluxes_at_phase).
        """
        i = int(np.argmin(np.abs(((self.rotation_phases - phase_deg + 180) % 360) - 180)))
        return Flux(self.wavelengths.copy(), self.fluxes[i, :].copy())

    def __str__(self):
        s = "Lightcurve Object:\n"
        s += f"  Phases (deg): N={self.rotation_phases.size}, range=({self.rotation_phases.min():.2f}, {self.rotation_phases.max():.2f})\n"
        s += f"  Wavelengths (μm): {self.wavelengths*1e6}\n"
        s += f"  Fluxes shape: {self.fluxes.shape}  [W/(m²·m)]\n"
        return s
