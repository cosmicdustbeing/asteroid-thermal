# shape_model_base.py
from abc import ABC, abstractmethod

class ShapeModelBase(ABC):
    @abstractmethod
    def get_facets(self):
        """
        Return a list of facet dictionaries.
        Each facet dictionary should include keys like:
            'area', 'normal', 'incidence_angle', 'emission_angle', etc.
        """
        pass

    
