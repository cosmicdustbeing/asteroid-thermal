# roastm.py
"""
This module implements the ROASTM thermal model as a class.
It accepts a ShapeModel instance and uses its facets’ solar incidence angles
to compute a temperature for each facet based on the ROASTM approach.
"""

import numpy as np
from flux import planck_radiance, Flux, Lightcurve
from constants import AU, S_SOLAR, SIGMA, H, C, K_B, PI

class ROASTM:
    def __init__(self, shape_model, A, f_iso, epsilon, r_au, eta_b = 1., ss_lon=0, ss_lat=0):
        """
        Initialize the ROASTM thermal model.
        
        Parameters:
            shape_model (ShapeModel): An instance containing facet data.
            A (float): Bond albedo (can be overall or per-facet).
            eta_b (float): Beaming parameter.
            epsilon (float): Emissivity (assumed constant for the model).
            r_au (float): Heliocentric distance in astronomical units.
            ss_lon (float, optional): Subsolar longitude in degrees (default: 0).
            ss_lat (float, optional): Subsolar latitude in degrees (default: 0).
        """
        self.shape_model = shape_model
        self.A = A
        self.eta_b = eta_b
        self.f_iso = f_iso
        self.epsilon = epsilon
        self.r_au = r_au
        self.r_m = r_au * AU
        self.ss_lon_rad = np.deg2rad(ss_lon)
        self.ss_lat_rad = np.deg2rad(ss_lat)

        # Compute the equilibrium, subsolar, and isolatitudinal temperatures.
        self.T_equ = self.equilibrium_temperature()

        self.T_ss = self.subsolar_temperature()
        
        self.T_iso = self.isothermal_temperature()

        # compute f_ss
        self.f_ss = self.fss_from_fiso()

        # maximum and minimum temperatures
        self.T_min = self.f_iso*self.T_iso
        self.T_max = self.T_min + self.f_ss*self.T_ss
        
        # time of T_max
        self.Tmax_rad = np.deg2rad(self.shift_Tmax())

        # Compute surface temperatures.
        self.T_surf = self.temperature_distribution()

    def fss_from_fiso(self):

        high = 0.516
        min0 = 0.997
        k = 0.134
        A = 1.000
        
        Tmax = A + (high - A) / (1 + np.exp((min0 - self.f_iso)/ k))

        f_ss = Tmax - self.f_iso*0.7552
        
        return f_ss

    def shift_Tmax(self):

        high = 47.05
        low = -0.02
        x_0 = 0.645
        k = 0.0868
        
        shift_lon = low + (high - low) / (1 + np.exp((x_0 - (self.f_iso*np.pi**-0.25)) / k))

        return shift_lon
    
    def equilibrium_temperature(self):
        """
        Compute the equilibrium temperature T_eq:

          T_eq = [ ( (1-A) * S_SOLAR ) / (epsilon * SIGMA * r_m^2) ]^(1/4)

        Return:
            float: Equilibrium subsolar temperature in Kelvin.        
        """
        numerator = (1 - self.A) * S_SOLAR
        denominator = self.epsilon * SIGMA * (self.r_au ** 2)

        T_eq = (numerator / denominator) ** 0.25
        return T_eq
        
    def subsolar_temperature(self):
        """
        Compute the subsolar temperature T_ss using the NEATM formula:

          T_ss = [ ( (1-A) * S_SOLAR ) / (eta * epsilon * SIGMA * r_m^2) ]^(1/4)

        Returns:
            float: Equilibrium subsolar temperature in Kelvin.
        """
        numerator = (1 - self.A) * S_SOLAR
        denominator = self.eta_b * self.epsilon * SIGMA * (self.r_au ** 2)

        T_ss = (numerator / denominator) ** 0.25
        return T_ss

    def isothermal_temperature(self):
        """
        Compute the latitudinal isothermal temperature T_iso using the FRM formula:

          T_ss_iso = [ ( (1-A) * S_SOLAR ) / (pi * epsilon * SIGMA * r_m^2) ]^(1/4)

        Returns:
            float: Isothermal subsolar temperature in Kelvin.
        """
        numerator = (1 - self.A) * S_SOLAR
        denominator = PI * self.epsilon * SIGMA * (self.r_au ** 2)

        T_iso = (numerator / denominator) ** 0.25
        return T_iso
    
    def temperature_distribution(self):
        """
        Compute the temperature for each facet in the shape model using
        spherical trigonometry. For each facet, the incidence angle is calculated
        using the facet's sun-normal latitude and longitude (sn_lat, sn_lon)
        and the subsolar coordinates (ss_lat, ss_lon) as follows:
        
            cos(theta) = sin(sn_lat) * sin(ss_lat) + cos(sn_lat) * cos(ss_lat) * cos(sn_lon - ss_lon)
        
        Then, if cos(theta) > 0, the facet temperature is assigned as:
        
            T = T_ss_max * (cos(theta))^(1/4)
        
        Otherwise, if cos(theta) <= 0 (i.e. facet is in shadow), T = 0.
        
        The computed temperatures are stored in the shape model's
        'temperatures' attribute.
        
        Returns:
        numpy array: Temperature for each facet.
        """

        facet_data = self.shape_model.facet_data
        sn_lat = facet_data['sn_lat']
        sn_lon = facet_data['sn_lon']

        cos_lat = np.sin(sn_lat) * np.sin(self.ss_lat_rad) + np.cos(sn_lat) * np.cos(self.ss_lat_rad)
        eff_cos_lat = np.where(cos_lat > 0, cos_lat, 0)
        
        cos_theta = np.sin(sn_lat) * np.sin(self.ss_lat_rad)
        cos_theta += np.cos(sn_lat) * np.cos(self.ss_lat_rad) * np.cos(sn_lon - self.ss_lon_rad - self.Tmax_rad)

        cos_theta = np.clip(cos_theta, -1.0, 1.0)
        eff_cos_theta = np.where(cos_theta > 0, cos_theta, 0)
        
        ss_temperatures = self.T_ss * (eff_cos_theta**0.25)
        iso_temperatures = self.T_iso * (eff_cos_lat**0.25)
        
        return (ss_temperatures * self.f_ss) + (iso_temperatures * self.f_iso)

    def roastm_flux(self, diam_eff, wavelength, Delta_au, so_lon=0, so_lat=0):
        """
        Compute the integrated flux from the ROASTM model at a given wavelength or array of wavelengths.
        Parameters:
            diam_eff (float): Effective diameter (in km) of the object.
            wavelength (float or numpy.ndarray): Wavelength(s) at which to compute the flux (in micrometers).
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
        cos_emission = np.sin(sn_lat) * np.sin(so_lat_rad)
        cos_emission += np.cos(sn_lat) * np.cos(so_lat_rad) * np.cos(sn_lon - so_lon_rad)

        # Eliminate cos_emission < 0.
        effective_cos = np.where(cos_emission > 0, cos_emission, 0)

        # Retrieve the temperatures computed by ROASTM.
        T = self.T_surf
        
        # wavelength is an array
        wavelengths = np.atleast_1d(wavelength)

        # Compute spectral radiance for each facet
        B_lambda = planck_radiance(wavelengths, T[:, None])

        # Compute the flux contribution from each facet:
        # facet_flux = epsilon * B_lambda * (facet_area * effective_cos)
        # Here, self.epsilon is used.
        facet_flux = self.epsilon * B_lambda * (facet_area * effective_cos)[:, None]

        # Sum flux contributions over all facets.
        F_total = np.sum(facet_flux, axis=0) * ((diam_eff / 2) ** 2) / (Delta_km ** 2)

        # Create a Flux instance to store the wavelengths and computed flux.
        flux_obj = Flux(wavelengths, F_total)

        return flux_obj

    def roastm_lightcurve(self, diam_eff, wavelength, Delta_au, n_steps, so_lon=0, so_lat=0):

        # Convert observer distance from AU to kilometers.
        Delta_km = Delta_au*AU*1e-3

        # Convert sub-observer coordinates to radians.
        so_lon_rad = np.deg2rad(so_lon)
        so_lat_rad = np.deg2rad(so_lat)

        # wavelength is an array
        wavelengths = np.atleast_1d(np.asarray(wavelength, dtype=float))
        n_lambda = wavelengths.size

        # calculate the rotation increment
        rotation_phases_deg = np.linspace(0.0, 360.0, int(n_steps), endpoint=False)
        d_rot_rad = np.deg2rad(360.0 / n_steps)  # rotate this much each step

        # initialize array for fluxes
        F_total_rot = np.zeros((n_steps, n_lambda), dtype=float)

        # for loop around the rotation phases
        for i in range(n_steps):

            # mirror roastm_flux for flux computation
            # Extract facet data from the shape model.
            facet_data = self.shape_model.facet_data
            sn_lat = facet_data['sn_lat']  # sun-normal latitudes (radians)
            sn_lon = facet_data['sn_lon']  # sun-normal longitudes (radians)
            facet_area = facet_data['area']  # facet areas

            # Compute cosine of the emission angle using spherical cosine law
            cos_emission = np.sin(sn_lat) * np.sin(so_lat_rad)
            cos_emission += np.cos(sn_lat) * np.cos(so_lat_rad) * np.cos(sn_lon - so_lon_rad)

            # Eliminate cos_emission < 0.
            effective_cos = np.where(cos_emission > 0, cos_emission, 0)

            # Retrieve the temperatures computed by ROASTM.
            T = self.T_surf

            # Compute spectral radiance for each facet
            B_lambda = planck_radiance(wavelengths, T[:, None])

            # Compute the flux contribution from each facet:
            # facet_flux = epsilon * B_lambda * (facet_area * effective_cos)
            # Here, self.epsilon is used.
            facet_flux = self.epsilon * B_lambda * (facet_area * effective_cos)[:, None]

            # Sum flux contributions over all facets at this rotation phase.
            F_total_rot[i, :] = np.sum(facet_flux, axis=0) * ((diam_eff / 2) ** 2) / (Delta_km ** 2)

            # rotate shape by changing the ss_lon and so_lon
            self.ss_lon_rad += d_rot_rad
            so_lon_rad += d_rot_rad
            self.T_surf = self.temperature_distribution()            

        # Create a Flux instance to store the wavelengths and computed flux.
        lightcurve_obj = Lightcurve(rotation_phases_deg, wavelengths, F_total_rot)

        return lightcurve_obj
            
        
    def __str__(self):
        s = "ROASTM Thermal Model:\n"
        s += f"  Equilibrium Temperature: {self.T_equ:.2f} K\n"
        s += f"  Subsolar Temperature: {self.T_ss:.2f} K\n"
        s += f"  Isothermal Temperature: {self.T_iso:.2f} K\n"
        s += f"  Maximum Temperature: {self.T_max:.2f} K\n"
        s += f"  Minimum Temperature: {self.T_min:.2f} K\n"
        s += f"  Sub-solar coordinates (lon, lat): {np.rad2deg(self.ss_lon_rad)}°, {np.rad2deg(self.ss_lat_rad)}°\n"
        s += f"  Bond albedo: {self.A}, Beaming parameter: {self.eta_b}, Emissivity: {self.epsilon}\n"
        s += f"  Heliocentric distance: {self.r_au} AU"
        return s
