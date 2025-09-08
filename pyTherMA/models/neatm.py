# neatm.py
"""
This module implements the NEATM thermal model as a class.
It accepts a ShapeModel instance and uses its facets’ solar incidence angles
to compute a temperature for each facet based on the NEATM prescription.
"""

import numpy as np
from flux import planck_radiance, Flux
from constants import AU, S_SOLAR, SIGMA, H, C, K_B, PI

class NEATM:
    def __init__(self, shape_model, A, eta, epsilon, r_au, ss_lon=0, ss_lat=0):
        """
        Initialize the NEATM thermal model.
        
        Parameters:
            shape_model (ShapeModel): An instance containing facet data.
            A (float): Bond albedo (can be overall or per-facet).
            eta (float): Beaming parameter.
            epsilon (float): Emissivity (assumed constant for the model).
            r_au (float): Heliocentric distance in astronomical units.
            ss_lon (float, optional): Subsolar longitude in degrees (default: 0).
            ss_lat (float, optional): Subsolar latitude in degrees (default: 0).
        """
        self.shape_model = shape_model
        self.A = A
        self.eta = eta
        self.epsilon = epsilon
        self.r_au = r_au
        self.r_m = r_au * AU
        self.ss_lon_rad = np.deg2rad(ss_lon)
        self.ss_lat_rad = np.deg2rad(ss_lat)
        
        # Compute the subsolar temperature.
        self.T_ss = self.subsolar_temperature()

        # Compute surface temperatures.
        self.T_surf = self.temperature_distribution()
        
    def subsolar_temperature(self):
        """
        Compute the subsolar temperature T_ss using the NEATM formula:
        
          T_ss = [ ( (1-A) * S_SOLAR ) / (eta * epsilon * SIGMA * r_m^2) ]^(1/4)
          
        Returns:
            float: Subsolar temperature in Kelvin.
        """
        numerator = (1 - self.A) * S_SOLAR
        denominator = self.eta * self.epsilon * SIGMA * (self.r_au ** 2)

        T_ss = (numerator / denominator) ** 0.25
        return T_ss

    def temperature_distribution(self):
        """
        Compute the temperature for each facet in the shape model using
        spherical trigonometry. For each facet, the incidence angle is calculated
        using the facet's sun-normal latitude and longitude (sn_lat, sn_lon)
        and the subsolar coordinates (ss_lat, ss_lon) as follows:
        
            cos(theta) = sin(sn_lat) * sin(ss_lat) + cos(sn_lat) * cos(ss_lat) * cos(sn_lon - ss_lon)
        
        Then, if cos(theta) > 0, the facet temperature is assigned as:
        
            T = T_ss * (cos(theta))^(1/4)
        
        Otherwise, if cos(theta) <= 0 (i.e. facet is in shadow), T = 0.
        
        The computed temperatures are stored in the shape model's
        'temperatures' attribute.
        
        Returns:
        numpy array: Temperature for each facet.
        """

        facet_data = self.shape_model.facet_data
        sn_lat = facet_data['sn_lat']
        sn_lon = facet_data['sn_lon']

        cos_theta = np.sin(sn_lat) * np.sin(self.ss_lat_rad)
        cos_theta += np.cos(sn_lat) * np.cos(self.ss_lat_rad) * np.cos(sn_lon - self.ss_lon_rad)

        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        effective_cos = np.where(cos_theta > 0, cos_theta, 0)
        
        temperatures = self.T_ss * (effective_cos**0.25)

        self.shape_model.temperatures = temperatures
        return temperatures


    def neatm_flux(self, diam_eff, wavelength, Delta_au, so_lon=0, so_lat=0):
        """
        Compute the integrated flux from the NEATM model at a given wavelength or array of wavelengths.
        Parameters:
            diam_eff (float): Effective diameter (in km) of the object.
            wavelength (float or numpy.ndarray): Wavelength(s) at which to compute the flux (in meters).
            Delta_au (float): Observer-centric distance in astronomical units.
            so_lon (float): Sub-observer longitude in degrees (default=0).
            so_lat (float): Sub-observer latitude in degrees (default=0).
        Returns:
            numpy.ndarray or float: Integrated flux (W/m²/m) for each wavelength.
        """

        # Convert observer distance from AU to kilometers.
        Delta_km = Delta_au*AU*1e-3

        # Convert sub-observer coordinates to radians.
        so_lon_rad = np.deg2rad(so_lon)
        so_lat_rad = np.deg2rad(so_lat)

        # Extract facet data from the shape model.
        facet_data = self.shape_model.facet_data
        sn_lat = facet_data['sn_lat']  # sun-normal latitudes (radians)
        sn_lon = facet_data['sn_lon']  # sun-normal longitudes (radians)
        facet_area = facet_data['area']  # facet areas
        
        # Compute cosine of the emission angle using spherical cosine law
        cos_emission = np.sin(sn_lat) * np.sin(so_lat)
        cos_emission += np.cos(sn_lat) * np.cos(so_lat) * np.cos(sn_lon - so_lon)

        # Eliminate cos_emission < 0.
        effective_cos = np.where(cos_emission > 0, cos_emission, 0)

        # Retrieve the temperatures computed by NEATM.
        T = self.T_surf
        
        # wavelength is an array
        wavelength = np.atleast_1d(wavelength)
        
        # Compute spectral radiance for each facet
        B_lambda = planck_radiance(wavelength, T[:, None])

        # Compute the flux contribution from each facet:
        # facet_flux = epsilon * B_lambda * (facet_area * effective_cos)
        # Here, self.epsilon is used.
        facet_flux = self.epsilon * B_lambda * (facet_area * effective_cos)[:, None]

        # Sum flux contributions over all facets.
        F_total = np.sum(facet_flux, axis=0) * ((diam_eff / 2) ** 2) / (Delta_km ** 2)

        # Create a Flux instance to store the wavelengths and computed flux.
        flux_obj = Flux(wavelength, F_total)

        return flux_obj

    def __str__(self):
        s = "NEATM Thermal Model:\n"
        s += f"  Subsolar Temperature: {self.T_ss:.2f} K\n"
        s += f"  Sub-solar coordinates (lon, lat): {np.rad2deg(self.ss_lon_rad)}°, {np.rad2deg(self.ss_lat_rad)}°\n"
        s += f"  Bond albedo: {self.A}, Beaming parameter: {self.eta}, Emissivity: {self.epsilon}\n"
        s += f"  Heliocentric distance: {self.r_au} AU"
        return s
