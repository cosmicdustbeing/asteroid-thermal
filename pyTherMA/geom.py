# calculate observing geometries (sub-solar latitude, longtidue) for each observation
# geometry.py
import numpy as np
from geometry_base import GeometryBase

def unit_vector(v):
    return v / np.linalg.norm(v)

class Geometry(GeometryBase):
    def __init__(self, shape, sun_xyz, observer_xyz, spin_params):
        """
        Initialize the geometry class.

        Parameters:
          shape: An instance of a shape model class that provides facet data via get_facets()
                 and (optionally) facet normals via get_normals(). (For example, an ellipsoid class.)
          sun_xyz: array-like (3,), the Sun position vector in the inertial or body-fixed frame.
          observer_xyz: array-like (3,), the observer position vector in the inertial or body-fixed frame.
          spin_params: a dictionary containing spin orientation parameters.
                       For example: {'spin_lat': ..., 'spin_lon': ...} (in radians)
                       (In this implementation, we assume that sun_xyz and observer_xyz are already
                        transformed to the body-fixed frame. If not, use spin_params to rotate them.)
        """
        self.shape = shape
        self.sun_xyz = np.array(sun_xyz, dtype=float)
        self.observer_xyz = np.array(observer_xyz, dtype=float)
        self.spin_params = spin_params  # May be used for future transformations

        # Current rotation phase (in radians) applied to the shape's body-centered longitudes.
        self.rotation_phase = 0.0

        # Compute initial geometry.
        self.compute_geometry()

    def compute_geometry(self):
        """
        Compute the incidence and emission angles for each facet.
        This function retrieves the facets from the shape and applies the current rotation phase to
        the body-centered longitudes before computing normals and angles.
        """
        # Retrieve the facet data (assumed structured array with 'body_lat' and 'body_lon')
        facets = self.shape.get_facets()
        n_facets = facets.shape[0]

        # Create a copy so as not to modify the original shape data.
        facets_rotated = facets.copy()
        # Update the body-centered longitude with the current rotation phase:
        facets_rotated['body_lon'] = (facets_rotated['body_lon'] + self.rotation_phase) % (2 * np.pi)
        
        # Try to get facet normals from the shape object.
        if hasattr(self.shape, 'get_normals'):
            # If normals are provided, assume they are for the unrotated shape and rotate them accordingly.
            # For a rotation about the z-axis, the rotation matrix is:
            c = np.cos(self.rotation_phase)
            s = np.sin(self.rotation_phase)
            Rz = np.array([[c, -s, 0],
                           [s,  c, 0],
                           [0,  0, 1]])
            normals = self.shape.get_normals()  # shape (n_facets, 3)
            # Apply rotation to the normals:
            normals = normals @ Rz.T
        else:
            # If normals are not provided, compute approximate normals from the rotated lat/lon.
            lat = facets_rotated['body_lat']
            lon = facets_rotated['body_lon']
            # For an ellipsoid, assume the shape object has attributes a, b, c.
            A = self.shape.a
            B = self.shape.b
            C = self.shape.c
            x = A * np.cos(lat) * np.cos(lon)
            y = B * np.cos(lat) * np.sin(lon)
            z = C * np.sin(lat)
            # The outward normal is proportional to (x/A^2, y/B^2, z/C^2)
            normals = np.vstack((x / (A**2), y / (B**2), z / (C**2))).T
            normals = normals / np.linalg.norm(normals, axis=1)[:, None]

        # Assume that sun_xyz and observer_xyz are already in the body-fixed frame.
        sun_dir = unit_vector(self.sun_xyz)
        obs_dir = unit_vector(self.observer_xyz)

        # Compute incidence angles: angle between facet normal and sun direction.
        dot_sun = np.einsum('ij,j->i', normals, sun_dir)
        dot_sun = np.clip(dot_sun, -1.0, 1.0)
        self.incidence_angles = np.arccos(dot_sun)  # in radians

        # Compute emission angles: angle between facet normal and observer direction.
        dot_obs = np.einsum('ij,j->i', normals, obs_dir)
        dot_obs = np.clip(dot_obs, -1.0, 1.0)
        self.emission_angles = np.arccos(dot_obs)  # in radians

        # Identify sub-solar and sub-observer facets as those with minimum incidence and emission angles.
        self.subsolar_facet = int(np.argmin(self.incidence_angles))
        self.subobserver_facet = int(np.argmin(self.emission_angles))

        # Optionally, store the rotated facets for later reference.
        self.rotated_facets = facets_rotated

    def rotate_shape(self, rotation_angle):
        """
        Rotate the shape by the given rotation_angle (in radians). This updates the internal rotation phase,
        and re-computes the geometry (i.e. re-calculates the sub-solar/sub-observer coordinates, incidence and emission angles).

        Parameters:
          rotation_angle: float, rotation angle in radians (added to the current rotation phase).
        """
        self.rotation_phase = (self.rotation_phase + rotation_angle) % (2 * np.pi)
        self.compute_geometry()

    def get_incidence_angles(self):
        return self.incidence_angles

    def get_emission_angles(self):
        return self.emission_angles

    def get_subsolar_facet(self):
        return self.subsolar_facet

    def get_subobserver_facet(self):
        return self.subobserver_facet
