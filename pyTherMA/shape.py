import numpy as np
from scipy.spatial import ConvexHull
from shape_model_base import ShapeModelBase

def fibonacci_sphere_points(n):
        """
        Generate n points on the unit sphere using the Fibonacci sphere algorithm.
        Parameters:
            n (int): Number of points to generate.
        Returns:
            numpy.ndarray: An (n x 3) array of points on the sphere.
        """
        indices = np.arange(0, n, dtype=float) + 0.5
        phi = np.arccos(1 - 2 * indices / n)
        golden_ratio = (1 + np.sqrt(5)) / 2
        theta = 2 * np.pi * indices / golden_ratio
        x = np.sin(phi) * np.cos(theta)
        y = np.sin(phi) * np.sin(theta)
        z = np.cos(phi)
        points = np.stack((x, y, z), axis=-1)
        return points


class Sphere(ShapeModelBase):
    def __init__(self, nlat, nlon):
        """
        Initialize a sphere shape with volume = 1 km^3.
        
        nlat, nlon: Number of grid bins for latitude and longitude.
        The sphere is discretized into facets on a grid in spherical coordinates.
        """
        # Compute the sphere radius from the volume.
        # Volume of a sphere: V = (4/3)*pi*R^3 = 1 km^3  =>  R = (3/(4*pi))^(1/3)
        self.radius = 1.
        self.volume = (4/3) * np.pi * self.radius**3  # should equal 1 km^3

        # Generate the grid of "sun-normal" coordinates.
        # For a sphere, we use latitude from -pi/2 to pi/2 and longitude from 0 to 2*pi.
        self.snlat = np.linspace(-np.pi/2, np.pi/2, nlat)
        self.snlon = np.linspace(0, 2*np.pi, nlon, endpoint=False)

        # Create a 2-D grid of these angles.
        snlon_img, snlat_img = np.meshgrid(self.snlon, self.snlat)

        # Store the sun-normal coordinates in a structured array.
        self.sncoords = np.rec.array(np.empty((nlat, nlon), 
                                     dtype=[('sn_lat', 'f8'), ('sn_lon', 'f8')]))
        self.sncoords.sn_lat = snlat_img
        self.sncoords.sn_lon = snlon_img

        # For a sphere centered at the origin, the body-centered coordinates
        # (the “true” geographic coordinates on the surface) are identical
        # to the sun-normal coordinates.
        self.latlon = self.sncoords

        # Compute the facet areas.
        # Here we approximate the sphere surface by dividing it into (nlat-1)*nlon quadrilateral facets.
        # For a sphere, the area of a spherical quadrilateral can be approximated as:
        #   area = R^2 * (sin(lat_edge_top) - sin(lat_edge_bottom)) * (lon_edge_width)
        lat_edges = np.linspace(-np.pi/2, np.pi/2, nlat)
        lon_edges = np.linspace(0, 2*np.pi, nlon+1)  # nlon+1 edges for nlon facets
        facet_areas = np.empty((nlat-1, nlon))
        for i in range(nlat-1):
            # Differential area in latitude is based on difference in sine of the edge latitudes.
            dS = self.radius**2 * (np.sin(lat_edges[i+1]) - np.sin(lat_edges[i]))
            dL = lon_edges[1] - lon_edges[0]  # uniform longitude spacing
            facet_areas[i, :] = np.abs(dS * dL)

        # Assemble the facet data into a structured array.
        # The number of facets is (nlat-1)*nlon.
        total_facets = (nlat - 1) * nlon
        dtype = np.dtype([
            ('body_lat', 'f8'),
            ('body_lon', 'f8'),
            ('sn_lat', 'f8'),
            ('sn_lon', 'f8'),
            ('area',   'f8')
        ])
        facets = np.empty(total_facets, dtype=dtype)
        facet_idx = 0
        # Compute facet centers from the edge coordinates.
        for i in range(nlat-1):
            center_lat = (lat_edges[i] + lat_edges[i+1]) / 2.0
            for j in range(nlon):
                center_lon = (lon_edges[j] + lon_edges[j+1]) / 2.0
                facets['body_lat'][facet_idx] = center_lat
                facets['body_lon'][facet_idx] = center_lon
                facets['sn_lat'][facet_idx] = center_lat  # same as body coordinates on a sphere
                facets['sn_lon'][facet_idx] = center_lon
                facets['area'][facet_idx] = facet_areas[i, j]
                facet_idx += 1

        self.facet_data = facets

        # -------------------------------
        # Compute facet vertices.
        # For each facet, the four vertices are determined by the edge coordinates.
        # We choose the following order:
        #   v1: (lat_edges[i+1], lon_edges[j])  -> top left
        #   v2: (lat_edges[i+1], lon_edges[j+1]) -> top right
        #   v3: (lat_edges[i],   lon_edges[j+1]) -> bottom right
        #   v4: (lat_edges[i],   lon_edges[j])   -> bottom left
        #
        # Convert spherical (r, lat, lon) to Cartesian (x, y, z):
        #   x = r * cos(lat) * cos(lon)
        #   y = r * cos(lat) * sin(lon)
        #   z = r * sin(lat)
        #
        facet_vertices = np.empty((total_facets, 4, 3), dtype=float)
        facet_idx = 0
        for i in range(nlat-1):
            for j in range(nlon):
                # Define the four vertices using the edge arrays.
                lat_bottom = lat_edges[i]
                lat_top = lat_edges[i+1]
                lon_left = lon_edges[j]
                lon_right = lon_edges[j+1]
                
                # Compute vertices in Cartesian coordinates.
                # Vertex order: top left, top right, bottom right, bottom left.
                v1 = self.radius * np.array([np.cos(lat_top) * np.cos(lon_left),
                                               np.cos(lat_top) * np.sin(lon_left),
                                               np.sin(lat_top)])
                v2 = self.radius * np.array([np.cos(lat_top) * np.cos(lon_right),
                                               np.cos(lat_top) * np.sin(lon_right),
                                               np.sin(lat_top)])
                v3 = self.radius * np.array([np.cos(lat_bottom) * np.cos(lon_right),
                                               np.cos(lat_bottom) * np.sin(lon_right),
                                               np.sin(lat_bottom)])
                v4 = self.radius * np.array([np.cos(lat_bottom) * np.cos(lon_left),
                                               np.cos(lat_bottom) * np.sin(lon_left),
                                               np.sin(lat_bottom)])
                facet_vertices[facet_idx, 0, :] = v1
                facet_vertices[facet_idx, 1, :] = v2
                facet_vertices[facet_idx, 2, :] = v3
                facet_vertices[facet_idx, 3, :] = v4
                facet_idx += 1

        self.facet_vertices = facet_vertices

        
    def get_facets(self):
        """Return the structured array of facets."""
        return self.facet_data

    def get_facet_vertices(self):
        """
        Returns the array of facet vertices.
        The returned array has shape (num_facets, 4, 3), where each facet is defined by four vertices
        in Cartesian coordinates.
        """
        return self.facet_vertices
    
    def get_volume(self):
        """Return the sphere volume (should equal 1 km^3)."""
        return self.volume

    def __str__(self):
        return f"SphereShape: Radius = {self.radius:.3f} km, {len(self.facet_data)} facets."

class FibonacciSphere(ShapeModelBase):
    def __init__(self, n_points=500):
        """
        Construct a sphere using the Fibonacci sphere method.
        This creates n_points approximately uniformly distributed over the unit sphere,
        and then computes the convex hull to generate triangular facets.

        Parameters:
            n_points (int): Number of points to generate on the sphere. The convex hull
                            will produce roughly 2*n_points - 4 facets (for n_points >= 4).
        """
        self.n_points = n_points

        # Volume of a sphere: V = (4/3)*pi*R^3 = 1 km^3  =>  R = (3/(4*pi))^(1/3)
        self.radius = 1.
        self.volume = (4/3) * np.pi * self.radius**3
        
        # Generate points on the sphere using the Fibonacci spiral.
        points = fibonacci_sphere_points(n_points)
        
        # Compute the convex hull to obtain the triangulation.
        hull = ConvexHull(points)
        n_facets = hull.simplices.shape[0]

            # Allocate arrays for facet vertices, areas, and centers.
        facet_vertices = np.empty((n_facets, 3, 3), dtype=float)
        facet_areas = np.empty(n_facets, dtype=float)
        facet_centers = np.empty((n_facets, 3), dtype=float)

        # Loop over each triangle facet from the convex hull.
        for i, simplex in enumerate(hull.simplices):
            v0, v1, v2 = points[simplex]
            facet_vertices[i, 0, :] = v0
            facet_vertices[i, 1, :] = v1
            facet_vertices[i, 2, :] = v2
            # Compute the area using the cross product (planar approximation).
            area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
            facet_areas[i] = area
            # Compute the center of the triangle.
            facet_centers[i] = (v0 + v1 + v2) / 3.0

        self.facet_vertices = facet_vertices

        # Convert each facet center from Cartesian to spherical coordinates (in radians).
        # For a point (x,y,z): lat = arcsin(z / r) and lon = arctan2(y, x)
        centers_r = np.linalg.norm(facet_centers, axis=1)
        body_lat = np.arcsin(facet_centers[:, 2] / centers_r)
        body_lon = np.arctan2(facet_centers[:, 1], facet_centers[:, 0])
        # Normalize longitudes to be in the range [0, 2pi).
        body_lon = np.mod(body_lon, 2 * np.pi)
        
        # For a sphere, the sun-normal coordinates are identical to the body-centered ones.
        sn_lat = body_lat.copy()
        sn_lon = body_lon.copy()

        # Assemble the facet data into a structured array.
        dtype = np.dtype([
            ('body_lat', 'f8'),
            ('body_lon', 'f8'),
            ('sn_lat', 'f8'),
            ('sn_lon', 'f8'),
            ('area', 'f8')
        ])
        facet_data = np.empty(n_facets, dtype=dtype)
        facet_data['body_lat'] = body_lat
        facet_data['body_lon'] = body_lon
        facet_data['sn_lat'] = sn_lat
        facet_data['sn_lon'] = sn_lon
        facet_data['area'] = facet_areas

        self.facet_data = facet_data

    def get_facet_vertices(self):
        """Return the array of facet vertices (shape: n_facets x 3 x 3)."""
        return self.facet_vertices

    def get_facets(self):
        """Return a structured array of facets."""
        return self.facet_data

    def get_volume(self):
        """Return the volume of the sphere (unit sphere volume)."""
        return self.volume

    def __str__(self):
        return f"FibonacciSphere: {self.faces.shape[0]} facets, Volume = {self.volume:.3f}"

    
class Ellipsoid(ShapeModelBase):
    def __init__(self, r1, r2, nlat, nlon):
        """
        r1: a/b axis ratio (user provided)
        r2: b/c axis ratio (user provided)
        nlat, nlon: number of bins for latitude and longitude
        """
        # Compute the semiaxes so that volume = 1 km^3.
        # Volume of ellipsoid = 4/3 * pi * a * b * c = 4/3 * pi * a^3 * (r2 / r1)
        # Set a^3 * (r2 / r1) = 1 => a = (r1/r2)^(1/3)
        self.a = (r1**2 * r2)**(1/3)
        self.b = self.a / r1
        self.c = self.b / r2  # or self.a/(r1*r2)

        # Compute the volume (for confirmation)
        self.volume = (4/3) * np.pi * self.a * self.b * self.c

        # Generate the grid of "sun-normal" coordinates analytically.
        # For example, we can use np.linspace to generate angles.
        # Note: The original IDL code computes snlat and snlon like:
        #    snlat = (findgen(nlat)/(nlat-1)*!dpi-!dpi/2.)
        #    snlon = ((findgen(nlon)*360/nlon+360.) mod 360.)*!dtor
        # In Python:
        self.snlat = np.linspace(-np.pi/2, np.pi/2, nlat)  # or appropriate formulation
        self.snlon = np.linspace(0, 2*np.pi, nlon, endpoint=False)

        # Create a meshgrid (or 2-D array) for these angles.
        snlon_img, snlat_img = np.meshgrid(self.snlon, self.snlat)

        # Compute the sun-normal coordinates (sncoords) as a structured 2D array.
        # In the IDL code, sncoords is a 2-d array of two rows: one for lat and one for lon.
        self.sncoords = np.rec.array(np.empty((nlat, nlon), dtype=[('sn_lat', 'f8'), ('sn_lon', 'f8')]))
        self.sncoords.sn_lat = snlat_img
        self.sncoords.sn_lon = snlon_img

        # Compute body-centered lat and lon (latlon) for each facet.
        # The IDL code uses transformations based on the ellipsoid axes.
        # For instance, it calculates:
        #     bctheta = atan(tan(snlon_img)*B/A)
        #     and then bcphi = atan(tan(snlat_img)*(C/A/B)*(Asinq+Bcosq)^(.5))
        # Here, A, B, C correspond to a, b, c.
        A, B, C = self.a, self.b, self.c

        # Compute bctheta:
        bctheta = np.arctan(np.tan(snlon_img) * (B/A))
        # There’s an adjustment in the IDL code to add pi offsets to some ranges.
        # (You’ll have to reproduce that logic; for brevity, we’ll note that a proper adjustment
        # is needed so that longitudes remain in [0, 2*pi] or [–pi, pi].)
        # Similarly, compute bcphi:
        Asinq = A**2 * np.sin(bctheta)**2
        Bcosq = B**2 * np.cos(bctheta)**2
        bcphi = np.arctan(np.tan(snlat_img) * (C/A/B) * np.sqrt(Asinq+Bcosq))

        # Create a structured array for body-centered coordinates:
        self.latlon = np.rec.array(np.empty((nlat, nlon), dtype=[('body_lat', 'f8'), ('body_lon', 'f8')]))
        self.latlon.body_lat = bcphi
        self.latlon.body_lon = bctheta

        # Compute the facet area for each facet.
        # The IDL code computes:
        #   dA = cos(bcphi)*sqrt( C^2*cos^2(bcphi)*(A^2*sin(theta)+B^2*cos(theta)) + A^2*B^2*sin^2(bcphi) )
        #   and then multiplies by dtheta and dphi (the grid differentials)
        # For our purposes, compute the grid differentials (dtheta, dphi) using np.diff and appropriate boundary fixes.
        # Then, compute the facet area as area = dA * dtheta * dphi.
        # You can follow a similar approach as in the IDL code.

        # (The details of the area computation can be placed in a separate method, e.g., self._compute_area())
        self.area = self._compute_area(snlat_img, snlon_img, bcphi, bctheta, nlat, nlon)

        # Store the facet data in a structured array.
        self.facet_data = self._assemble_facet_data()

    def _compute_area(self, snlat_img, snlon_img, bcphi, bctheta, nlat, nlon):
        # Compute the differential increments in theta (longitude) and phi (latitude)
        # For simplicity, assume uniform spacing:
        dtheta = (2*np.pi) / nlon
        dphi = np.pi / (nlat-1)
        # Following the IDL logic, the differential area element for a facet on an ellipsoid:
        #    dA = cos(bcphi) * sqrt( C^2*cos^2(bcphi)*(A^2*sin(bctheta) + B^2*cos(bctheta)) + A^2*B^2*sin^2(bcphi) )
        # Note: You’ll need to adjust this formula to match the IDL exactly.
        A, B, C = self.a, self.b, self.c
        # We create a function that computes the local area element:
        dA = np.cos(bcphi) * np.sqrt(
            C**2 * np.cos(bcphi)**2 * (A**2 * np.sin(bctheta)**2 + B**2 * np.cos(bctheta)**2) +
            (A * B)**2 * np.sin(bcphi)**2
        )
        # Multiply by dtheta and dphi to get the facet area:
        area = dA * dtheta * dphi
        return area

    def _assemble_facet_data(self):
        """
        Create a structured numpy array with fields:
        'body_lat', 'body_lon', 'sn_lat', 'sn_lon', 'area'
        The arrays (for each facet) should be flattened to 1D.
        """
        nlat, nlon = self.sncoords.sn_lat.shape
        total_facets = nlat * nlon
        dtype = np.dtype([
            ('body_lat', 'f8'),
            ('body_lon', 'f8'),
            ('sn_lat', 'f8'),
            ('sn_lon', 'f8'),
            ('area', 'f8')
        ])
        facets = np.empty(total_facets, dtype=dtype)
        facets['body_lat'] = self.latlon.body_lat.flatten()
        facets['body_lon'] = self.latlon.body_lon.flatten()
        facets['sn_lat'] = self.sncoords.sn_lat.flatten()
        facets['sn_lon'] = self.sncoords.sn_lon.flatten()
        facets['area'] = self.area.flatten()
        return facets

    def get_facets(self):
        """Return the structured array of facets."""
        return self.facet_data

    def get_volume(self):
        """Return the ellipsoid volume (should be 1 km^3)."""
        return self.volume

    def __str__(self):
        s = f"EllipsoidShape: semi-axes: a = {self.a:.3f},\n"
        s += f"a = {self.b:.3f}, c = {self.c:.3f}/n"
        s += f"{len(self.facet_data)} facets."
        return s

class FibonacciEllipsoid(ShapeModelBase):
    def __init__(self, r1, r2, n_points):
        """
        Initialize an ellipsoid shape using the Fibonacci algorithm.
        
        Parameters:
            r1 (float): a/b axis ratio (user provided).
            r2 (float): b/c axis ratio (user provided).
            n_points (int): Number of points to generate using the Fibonacci sphere algorithm.
            
        The semiaxes are computed so that the ellipsoid volume is 1 km^3.
        """
        # Compute semiaxes so that the volume is 1 km^3.
        # Using: volume = (4/3)*pi * a * b * c and
        #           a = (r1^2 * r2)^(1/3),  b = a / r1,  c = b / r2  (i.e. a/(r1*r2))
        self.a = (r1**2 * r2)**(1/3)
        self.b = self.a / r1
        self.c = self.b / r2  # equivalently, self.a/(r1*r2)
        
        # Compute the ellipsoid volume (for confirmation)
        self.volume = (4/3) * np.pi * self.a * self.b * self.c

        # Generate points on a unit sphere using the Fibonacci algorithm.
        points_unit = fibonacci_sphere_points(n_points)
        # Scale the points to form an ellipsoid: x' = a*x, y' = b*y, z' = c*z.
        ellipsoid_points = np.empty_like(points_unit)
        ellipsoid_points[:, 0] = self.a * points_unit[:, 0]
        ellipsoid_points[:, 1] = self.b * points_unit[:, 1]
        ellipsoid_points[:, 2] = self.c * points_unit[:, 2]
        self.points = ellipsoid_points

        # Compute the convex hull of the scaled points to obtain triangular facets.
        hull = ConvexHull(ellipsoid_points)
        n_facets = hull.simplices.shape[0]

        # Allocate arrays for facet vertices, centers, and areas.
        facet_vertices = np.empty((n_facets, 3, 3), dtype=float)
        facet_centers = np.empty((n_facets, 3), dtype=float)
        facet_areas = np.empty(n_facets, dtype=float)

        for i, simplex in enumerate(hull.simplices):
            v0, v1, v2 = ellipsoid_points[simplex]
            facet_vertices[i, 0, :] = v0
            facet_vertices[i, 1, :] = v1
            facet_vertices[i, 2, :] = v2
            # Calculate triangle area via cross product (planar approximation)
            area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
            facet_areas[i] = area
            # Compute the triangle (facet) center as the average of its vertices.
            facet_centers[i] = (v0 + v1 + v2) / 3.0

        self.facet_vertices = facet_vertices

        # Convert each facet center to spherical coordinates.
        # For a point (x, y, z):  latitude = arcsin(z / r_eff) and longitude = arctan2(y, x)
        centers_r = np.linalg.norm(facet_centers, axis=1)
        body_lat = np.arcsin(facet_centers[:, 2] / centers_r)
        body_lon = np.arctan2(facet_centers[:, 1], facet_centers[:, 0])
        # Normalize longitudes to the range [0, 2*pi)
        body_lon = np.mod(body_lon, 2 * np.pi)

        # In this implementation the sun-normal coordinates are taken to be identical.
        sn_lat = body_lat.copy()
        sn_lon = body_lon.copy()

        # Assemble the facet data into a structured NumPy array.
        dtype = np.dtype([
            ('body_lat', 'f8'),
            ('body_lon', 'f8'),
            ('sn_lat', 'f8'),
            ('sn_lon', 'f8'),
            ('area',   'f8')
        ])
        facet_data = np.empty(n_facets, dtype=dtype)
        facet_data['body_lat'] = body_lat
        facet_data['body_lon'] = body_lon
        facet_data['sn_lat'] = sn_lat
        facet_data['sn_lon'] = sn_lon
        facet_data['area'] = facet_areas

        self.facet_data = facet_data

        # For this Fibonacci-based shape, the full "grid" of sun-normal coordinates is not available.
        # Instead, we store the per-facet center coordinates in facet_data.
        self.sncoords = None
        self.latlon = None

    def get_facets(self):
        """Return the structured array of facet data."""
        return self.facet_data

    def get_facet_vertices(self):
        """Return the array of facet vertices (shape: [n_facets, 3, 3])."""
        return self.facet_vertices

    def get_volume(self):
        """Return the ellipsoid volume."""
        return self.volume

    def __str__(self):
        return (f"FibonacciEllipsoid: semiaxes: a = {self.a:.3f}, b = {self.b:.3f}, c = {self.c:.3f}, "
                f"{len(self.facet_data)} facets.")


class FibonacciCellinoid(ShapeModelBase):
    def __init__(self, r1, r2, ra, rb, rc, n_points):
        """
        Initialize an ellipsoid shape using the Fibonacci algorithm.
        
        Parameters:
            r1 (float): (a1+a2)/(b1+b2) axis ratio (user provided).
            r2 (float): (b1+b2)/(c1+c2) axis ratio (user provided).
            ra (float): a2/a1 axis ratio.
            rb (float): b2/b1 axis ratio.
            rc (float): c2/c1 axis ratio.
            n_points (int): Number of points to generate using the Fibonacci sphere algorithm.

        The semiaxes are computed so that the ellipsoid volume is 1 km^3.
        """
        # Compute semiaxes so that the volume is 1 km^3.
        # Using: volume = (4/3)*pi * a * b * c and
        #           a = (r1^2 * r2)^(1/3),  b = a / r1,  c = b / r2  (i.e. a/(r1*r2))
        self.a = (r1**2 * r2)**(1/3)
        self.b = self.a / r1
        self.c = self.b / r2  # equivalently, self.a/(r1*r2)
        
        # Compute the ellipsoid volume (for confirmation)
        self.volume = (4/3) * np.pi * self.a * self.b * self.c

        # Compute the asymmetric scaling factors:
        # For the x-axis:
        self.a1 = (2 * self.a) / (1 + ra)   # scaling factor for x >= 0
        self.a2 = (2 * self.a * ra) / (1 + ra)  # scaling factor for x < 0
        
        # For the y-axis:
        self.b1 = (2 * self.b) / (1 + rb)
        self.b2 = (2 * self.b * rb) / (1 + rb)
        
        # For the z-axis:
        self.c1 = (2 * self.c) / (1 + rc)
        self.c2 = (2 * self.c * rc) / (1 + rc)
        
        # Generate points on a unit sphere using the Fibonacci algorithm.
        points_unit = fibonacci_sphere_points(n_points)
        # Scale the points to form an ellipsoid: x' = a*x, y' = b*y, z' = c*z.
        ellipsoid_points = np.empty_like(points_unit)

        # For x-coordinate: if x is greater than or equal to zero, scale by a1; otherwise, scale by a2.
        ellipsoid_points[:, 0] = np.where(points_unit[:, 0] >= 0,
                                          self.a1 * points_unit[:, 0],
                                          self.a2 * points_unit[:, 0])

        # For y-coordinate: if y is greater than or equal to zero, scale by b1; otherwise, scale by b2.
        ellipsoid_points[:, 1] = np.where(points_unit[:, 1] >= 0,
                                          self.b1 * points_unit[:, 1],
                                          self.b2 * points_unit[:, 1])

        # For z-coordinate: if z is greater than or equal to zero, scale by c1; otherwise, scale by c2.
        ellipsoid_points[:, 2] = np.where(points_unit[:, 2] >= 0,
                                          self.c1 * points_unit[:, 2],
                                          self.c2 * points_unit[:, 2])
        
        self.points = ellipsoid_points

        # Compute the convex hull of the scaled points to obtain triangular facets.
        hull = ConvexHull(ellipsoid_points)
        n_facets = hull.simplices.shape[0]

        # Allocate arrays for facet vertices, centers, and areas.
        facet_vertices = np.empty((n_facets, 3, 3), dtype=float)
        facet_centers = np.empty((n_facets, 3), dtype=float)
        facet_areas = np.empty(n_facets, dtype=float)

        for i, simplex in enumerate(hull.simplices):
            v0, v1, v2 = ellipsoid_points[simplex]
            facet_vertices[i, 0, :] = v0
            facet_vertices[i, 1, :] = v1
            facet_vertices[i, 2, :] = v2
            # Calculate triangle area via cross product (planar approximation)
            area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
            facet_areas[i] = area
            # Compute the triangle (facet) center as the average of its vertices.
            facet_centers[i] = (v0 + v1 + v2) / 3.0

        self.facet_vertices = facet_vertices

        # Convert each facet center to spherical coordinates.
        # For a point (x, y, z):  latitude = arcsin(z / r_eff) and longitude = arctan2(y, x)
        centers_r = np.linalg.norm(facet_centers, axis=1)
        body_lat = np.arcsin(facet_centers[:, 2] / centers_r)
        body_lon = np.arctan2(facet_centers[:, 1], facet_centers[:, 0])
        # Normalize longitudes to the range [0, 2*pi)
        body_lon = np.mod(body_lon, 2 * np.pi)

        # In this implementation the sun-normal coordinates are taken to be identical.
        sn_lat = body_lat.copy()
        sn_lon = body_lon.copy()

        # Assemble the facet data into a structured NumPy array.
        dtype = np.dtype([
            ('body_lat', 'f8'),
            ('body_lon', 'f8'),
            ('sn_lat', 'f8'),
            ('sn_lon', 'f8'),
            ('area',   'f8')
        ])
        facet_data = np.empty(n_facets, dtype=dtype)
        facet_data['body_lat'] = body_lat
        facet_data['body_lon'] = body_lon
        facet_data['sn_lat'] = sn_lat
        facet_data['sn_lon'] = sn_lon
        facet_data['area'] = facet_areas

        self.facet_data = facet_data

        # For this Fibonacci-based shape, the full "grid" of sun-normal coordinates is not available.
        # Instead, we store the per-facet center coordinates in facet_data.
        self.sncoords = None
        self.latlon = None

    def get_facets(self):
        """Return the structured array of facet data."""
        return self.facet_data

    def get_facet_vertices(self):
        """Return the array of facet vertices (shape: [n_facets, 3, 3])."""
        return self.facet_vertices

    def get_volume(self):
        """Return the ellipsoid volume."""
        return self.volume

    def __str__(self):
        return (f"FibonacciEllipsoid: semiaxes: a = {self.a:.3f}, b = {self.b:.3f}, c = {self.c:.3f}, "
                f"{len(self.facet_data)} facets.")
