# plotting.py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.colors as colors

def plot_shape_mesh(shape, 
                    title="Shape Mesh",
                    edge_color='k',
                    face_color='grey',
                    alpha=0.4):
    """
    Plot a 3D mesh of a shape's facets using facet vertices.
    
    The shape object must have an attribute `facet_vertices` that is a NumPy array 
    of shape (n_facets, n_vertices, 3), where each facet is defined by the (x, y, z)
    coordinates of its vertices.
    
    This function works for any shape model including those loaded from .obj files.
    
    Parameters:
        shape: An instance of a shape class with a 'facet_vertices' attribute.
        title (str): The title for the plot.
        edge_color (str): The color for the facet edges.
        face_color (str): The color for the facet faces.
        alpha (float): The transparency of the facet faces (0 is transparent, 1 is opaque).
    """
    # Check if shape has facet_vertices
    if not hasattr(shape, 'facet_vertices'):
        raise AttributeError("The provided shape does not have a 'facet_vertices' attribute.")

    facets = shape.facet_vertices  # Expected shape: (n_facets, n_vertices, 3)
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal', 'box')
    
    # Create a Poly3DCollection from the facet vertices.
    poly_collection = Poly3DCollection(facets, 
                                       facecolors=face_color, 
                                       edgecolors=edge_color, 
                                       linewidths=0.2,
                                       alpha=alpha)
    ax.add_collection3d(poly_collection)
    
    # Set plot limits by finding the extents of the vertices.
    all_vertices = facets.reshape(-1, 3)
    x_min, y_min, z_min = np.min(all_vertices, axis=0)
    x_max, y_max, z_max = np.max(all_vertices, axis=0)
    
    # Expand limits slightly for better visualization.
    padding = 0.1 * (x_max - x_min)
    ax.set_xlim([x_min - padding, x_max + padding])
    ax.set_ylim([y_min - padding, y_max + padding])
    ax.set_zlim([z_min - padding, z_max + padding])
    
    ax.set_title(title)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    
    plt.show()

def plot_thermal_mesh(shape, thermal_model, 
                      title="Surface Temperature Map", 
                      cmap_name="rainbow",
                      vmin=None, vmax=None):
    """
    Plot a 3D mesh of the shape facets colored by their temperatures.
    
    Parameters:
        shape: A shape model instance that has an attribute `facet_vertices`
               (an array of shape (n_facets, n_vertices, 3)) containing the vertices
               of each facet.
        thermal_model: A thermal model instance that contains an attribute
                       T_surf, an array of temperatures (length = n_facets) for each facet.
        title (str): Title for the plot.
        cmap_name (str): Name of the Matplotlib colormap to use.
        vmin (float): Minimum temperature for colormap scaling. If None, uses T.min().
        vmax (float): Maximum temperature for colormap scaling. If None, uses T.max().
    
    Returns:
        None
    """
    # Retrieve facet vertices from the shape.
    if not hasattr(shape, "facet_vertices"):
        raise AttributeError("The provided shape does not have a 'facet_vertices' attribute.")
    facets = shape.facet_vertices  # shape: (n_facets, n_vertices, 3)
    
    # Retrieve temperatures from the thermal model.
    T = np.array(thermal_model.T_surf)  # shape: (n_facets,)
    
    # If vmin or vmax are not provided, set them based on the temperature array.
    if vmin is None:
        vmin = T.min()
    if vmax is None:
        vmax = T.max()
    
    # Normalize temperatures for the colormap.
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)
    facecolors = cmap(norm(T))
    
    # Create a Poly3DCollection with facet facecolors.
    poly_collection = Poly3DCollection(facets, facecolors=facecolors, edgecolors='k', alpha=0.9)
    
    # Set up the 3D plot.
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_aspect('equal', 'box')
    ax.add_collection3d(poly_collection)
    
    # Set plot limits using all vertices.
    all_vertices = facets.reshape(-1, 3)
    x_min, y_min, z_min = np.min(all_vertices, axis=0)
    x_max, y_max, z_max = np.max(all_vertices, axis=0)
    padding = 0.1 * (x_max - x_min)
    ax.set_xlim([x_min - padding, x_max + padding])
    ax.set_ylim([y_min - padding, y_max + padding])
    ax.set_zlim([z_min - padding, z_max + padding])
    
    ax.set_title(title)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    
    # Add a colorbar. We create a mappable object for the colorbar.
    mappable = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array(T)
    cbar = plt.colorbar(mappable, ax=ax, pad=0.1)
    cbar.set_label("Temperature (K)")
    
    plt.show()

def plot_temperature_map(shape, thermal_model,
                         title="Temperature Contour Map", cmap="rainbow"):
    """
    Plot a 2D contour map of the surface temperatures from the shape and thermal model.
    
    Parameters:
        shape: A shape model instance with an attribute `facet_data` containing at least
               'body_lat' and 'body_lon' (in radians) for each facet.
        thermal_model: A thermal model instance that contains an attribute (e.g., T_surf)
                       which is an array of temperatures corresponding to each facet.
        title (str): Title of the plot.
        cmap (str): Matplotlib colormap to use.
    """
    # Extract facet latitudes and longitudes (in radians) from the shape model.
    facets = shape.facet_data
    lat = np.array(facets['body_lat'])
    lon = np.array(facets['body_lon'])
    
    # Extract the temperature for each facet.
    T = np.array(thermal_model.T_surf)
    
    # Convert lat and lon from radians to degrees for plotting.
    lat_deg = np.degrees(lat)
    lon_deg = np.degrees(lon)
    
    # Create a tricontourf plot for scattered data.
    plt.figure(figsize=(8, 6))
    contour = plt.tricontourf(lon_deg, lat_deg, T, levels=20, cmap=cmap)
    plt.scatter(np.rad2deg(thermal_model.ss_lon_rad),np.rad2deg(thermal_model.ss_lat_rad),
                marker='X',color='black', label='Sub-Solar Point')
    plt.colorbar(contour, label="Temperature (K)")
    plt.title(title)
    plt.xlabel("Longitude (deg)")
    plt.ylabel("Latitude (deg)")
    plt.tight_layout()
    plt.legend()
    plt.show()

def plot_lightcurve(lc, ax=None, unit="W/(m²·m)", to_Jy=False, **kwargs):
    """
    Plot the lightcurve flux vs rotation phase.

    Parameters
    ----------
    lc : Lightcurve
        The Lightcurve instance to plot.
    ax : matplotlib.axes.Axes, optional
        Axis to plot on. If None, a new figure+axes is created.
    unit : str
        Label for the y-axis if plotting in SI units.
    to_Jy : bool
        If True, convert to Jansky.
    **kwargs :
        Extra keyword arguments passed to matplotlib plot().
    """
    if ax is None:
        fig, ax = plt.subplots()

    # Select flux array
    if to_Jy:
        flux = lc.to_Jy()
        y_label = "Flux density (Jy)"
    else:
        flux = lc.fluxes
        y_label = f"Flux density [{unit}]"

    # Loop wavelengths
    for j, lam in enumerate(lc.wavelengths):
        lam_um = lam  # convert m → µm
        ax.plot(
            lc.rotation_phases,
            flux[:, j],
            label=f"{lam_um:.2f} µm",
            **kwargs
        )

    ax.set_xlabel("Rotation phase (deg)")
    ax.set_ylabel(y_label)
    ax.legend()
    ax.set_xlim(0, 360)
    plt.show()
