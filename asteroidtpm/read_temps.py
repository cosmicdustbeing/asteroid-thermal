#!/usr/bin/env /usr/local/opt/python@3.10/bin/python3.10

import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

A_b = 0.03  # Bond albedo
R_au = 2.0 # heliocentric distance in AU
emis = 0.9

# sub-solar equilibrium temperature
tss = (1367.*(1-A_b)/(emis*5.67e-8*R_au*R_au))**.25

# load smooth surface temperatures
data = np.load('surftemp_smooth_sslat60_theta2.npz')

# set up smooth surface coordinates
lat=np.linspace(-90., 90., data['smooth'].shape[0],dtype=float)
lon=np.arange(0.,360.,360./data['smooth'].shape[1],dtype=float)

# data array containing smooth temperatures
da = xr.DataArray(data['smooth'], coords=[lat, lon], dims=["latitude", "longitude"])

# data set
ds = xr.Dataset({"smooth":da})

# load crater temperatures
data = np.load('surftemp_crater68_sslat60_theta2.npz')

# set up cratered surface coordinates
ds.coords["crlat"] = np.linspace(-90., 90., data['crel1'].shape[0],dtype=float)
ds.coords["crlon"] = np.arange(0.,360.,360./data['crel1'].shape[1],dtype=float)

# loop through crater elements and store in ds data structure
for key in data:
    if key != "smooth":
        ds[key]=(['crlat', 'crlon'],  data[key])

# scale array by sub-solar temperature
ds = ds*tss

print(ds.max())


# interpolate across latitudes and plot smooth surface tepmeratures
sm_interp = ds['smooth'].interp(latitude=np.linspace(-90., 90.,100,dtype=float),method="quadratic")
sm_interp.where(sm_interp > 10., other = 10.).plot.contourf(cmap=mpl.cm.RdYlBu_r)

#ds['smooth'].interp(latitude=np.linspace(-90., 90.,100,dtype=float),method="quadratic").plot()
#ds['crel4'].interp(crlat=np.linspace(-90., 90.,100,dtype=float)).plot()
plt.show()
