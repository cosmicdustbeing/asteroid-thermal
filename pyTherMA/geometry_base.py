# geometry_base.py
from abc import ABC, abstractmethod

class GeometryBase(ABC):
    @abstractmethod
    def compute_geometry(self):
        """
        Compute the geometry for the current epoch.
        This should compute, for each facet:
          - the incidence angle (angle between facet normal and sun direction)
          - the emission angle (angle between facet normal and observer direction)
          - (optionally) identify the sub-solar and sub-observer facets.
        """
        pass

    @abstractmethod
    def get_incidence_angles(self):
        """Return the incidence angles (in radians) as a 1-D numpy array."""
        pass

    @abstractmethod
    def get_emission_angles(self):
        """Return the emission angles (in radians) as a 1-D numpy array."""
        pass

    @abstractmethod
    def get_subsolar_facet(self):
        """Return the index (or indices) of the facet(s) closest to the sub-solar point."""
        pass

    @abstractmethod
    def get_subobserver_facet(self):
        """Return the index (or indices) of the facet(s) closest to the sub-observer point."""
        pass

    @abstractmethod
    def rotate_shape(self, rotation_angle):
        """
        Rotate the shape by the given rotation_angle (in radians) in the longitude direction.
        This updates the geometry by recomputing the sub-solar/sub-observer coordinates,
        incidence and emission angles, etc.
        """
        pass
