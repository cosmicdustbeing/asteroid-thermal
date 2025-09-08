# iso.py
"""
This module implements the Isothermal thermal model as a class.
It accepts a ShapeModel instance and uses its facets’ solar incidence angles
to compute a temperature for each facet based on the isothermal model.
"""

import numpy as np
from flux import planck_radiance, Flux
from constants import AU, S_SOLAR, SIGMA, H, C, K_B, PI

class ISO:
    def __init__(self, shape_model, A, epsilon, r_au):
        """
        Initialize the Isothermal model.
        
        Parameters:
            shape_model (ShapeModel): An instance containing facet data.
            A (float): Bond albedo (can be overall or per-facet).
            epsilon (float): Emissivity (assumed constant for the model).
            r_au (float): Heliocentric distance in astronomical units.
        """
        self.shape_model = shape_model
        self.A = A
        self.epsilon = epsilon
        self.r_au = r_au
        self.r_m = r_au * AU

        # Compute the isothermal temperature.
        self.T_iso = self.isothermal_temperature()

        self.T_surf = self.temperature_distribution()
        
    def isothermal_temperature(self):
        """
        Compute the subsolar temperature T_iso using the isothermal:
        T_iso = [ ( (1-A) * S_SOLAR ) / (4 * pi * epsilon * SIGMA * r_m^2) ]^(1/4)          
        Returns:
        float: Isothermal temperature in Kelvin.
        """
        numerator = (1 - self.A) * S_SOLAR
        denominator = 4 * PI * self.epsilon * SIGMA * (self.r_au ** 2)

        T_iso = (numerator / denominator) ** 0.25
        return T_iso

    def temperature_distribution(self):
        """
        Compute the temperature for each facet in the shape model.
        For each facet (assumed to have solar-relative latitude and longitude
        stored as 'lat' and 'lon' in the facet dictionary), we assign:
        
        T = T_iso for all facets,
        The computed temperatures are stored in the shape model's 'temperatures' attribute.
    
        Returns:
        numpy array: Temperature for each facet.
        """
        
        n_facets = len(self.shape_model.facet_data)

        # Create a flat array filled with T_iso
        temperatures = np.full(n_facets, self.T_iso, dtype=float)

        # Store on the shape model for downstream use
        self.shape_model.temperatures = temperatures

        return temperatures

    def iso_flux(self, diam_eff, wavelength, Delta_au, so_lat=0):
        """
        Compute the integrated flux from the Isothermal model at a given wavelength or array of wavelengths.
        Parameters:
            diam_eff (float): Effective diameter (in km) of the object.
            wavelength (float or numpy.ndarray): Wavelength(s) at which to compute the flux (in meters).
            Delta_au (float): Observer-centric distance in astronomical units.
            so_lat (float): Sub-observer latitude in degrees (default=0).
        Returns:
            numpy.ndarray or float: Integrated flux (W/m²/m) for each wavelength.
        """
        
        # Convert observer distance from AU to kilometers.
        Delta_km = Delta_au*AU*1e-3

        # Convert sub-observer coordinates to radians.
        so_lat_rad = np.deg2rad(so_lat)

        # Extract facet data from the shape model.
        facet_data = self.shape_model.facet_data
        sn_lat = facet_data['sn_lat']  # sun-normal latitudes (radians)
        sn_lon = facet_data['sn_lon']  # sun-normal longitudes (radians)
        facet_area = facet_data['area']  # facet areas

        # Compute cosine of the emission angle using spherical cosine law
        cos_emission = np.sin(sn_lat) * np.sin(so_lat)
        cos_emission += np.cos(sn_lat) * np.cos(so_lat) * np.cos(sn_lon)

        # Eliminate cos_emission < 0.
        effective_cos = np.where(cos_emission > 0, cos_emission, 0)

        # Retrieve the temperatures computed by ISO.
        T = self.T_surf

        # wavelength is an array
        wavelength = np.atleast_1d(wavelength)

        # Compute spectral radiance for each facet
        B_lambda = planck_radiance(wavelength, T[:, None])

        # Compute the flux contribution from each facet:
        facet_flux = self.epsilon * B_lambda * (facet_area * effective_cos)[:, None]

        # Sum flux contributions over all facets.
        F_total = np.sum(facet_flux, axis=0) * ((diam_eff / 2) ** 2) / (Delta_km ** 2)

        # Create a Flux instance to store the wavelengths and computed flux.
        flux_obj = Flux(wavelength, F_total)

        return flux_obj
        
    def __str__(self):
        s = "ISO Thermal Model:\n"
        s += f"  Isothermal Temperature: {self.T_iso:.2f} K\n"
        s += f"  Bond albedo: {self.A}, Emissivity: {self.epsilon}\n"
        s += f"  Heliocentric distance: {self.r_au} AU"
        return s
