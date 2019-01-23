from __future__ import division
import sys
sys.path.append('/usr/projects/climate/jrub/pyfunctions')
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, addcyclic
from funciones import *
from classes import *
matplotlib.use("Agg")
import csv
import xarray as xr

#-------
# X = int(sys.argv[1]) #200
# nlayer = int(sys.argv[2]) #300
# ne = int(sys.argv[3]) #60
# hydro = int(sys.argv[4]) #0
# lon0 = int(sys.argv[5]) #270
# dframe = int(sys.argv[6]) #10
# tfile_init = sys.argv[7]
# tfile_fin  = sys.argv[8]
# ------
X = 200
nlayer = 30
ne = 60
hydro = 1
lon0 = 0 #0
dframe = 1
tfile_init = 0
tfile_fin  = 17
#-------


variables = ('zeta_i','gradTh_i','Va_gTh_i','zeta_j','gradTh_j','Va_gTh_j','zeta_k','gradTh_k','Va_gTh_k','u','Th','Va_gTh')

colmaps = ('seismic','seismic','seismic','seismic','seismic','seismic','seismic','seismic','seismic',
           'seismic','jet','seismic')


label =(r"$\xi_{i}$",r"$\Theta_{i}$",r"$\xi_{i} \cdot \Theta_{i}$",
        r"$\xi_{j}$",r"$\Theta_{j}$",r"$\xi_{j} \cdot \Theta_{j}$", 
        r"$\xi_{k}$",r"$\Theta_{k}$",r"$(\xi_{k}+f)\cdot \Theta_{k}$", 
        r"u",r"$\Theta$",r"$(\xi + f) \cdot \nabla \Theta$")

# limits =([-3e-9,3e-9],
#          [-1e-4,1e-4],
#          [-5e-2,5e-2],

#          [-3e-5,3e-5],                                                                                     
#          [-0.007,0.007],
#          [-5e-2,5e-2],

#          [-3e-3,3e-3],
#          [-8e-2,8e-2],
#          [-5e-3,5e-3],

#          [-30,30],
#          [240,820],
#          [-5e-2,5e-2])

limits =([-3e-9,3e-9],
         [-1e-3,1e-3],
         [-5e-2,5e-2],

         [-3e-2,3e-2],                                                                                     
         [-0.007,0.007],
         [-5e-2,5e-2],

         [-3e-3,3e-3],
         [-8e-2,8e-2],
         [-5e-3,5e-3],

         [-30,30],
         [240,820],
         [-5e-2,5e-2])

if hydro :
    dyn = 'hydro'
else:
    dyn = 'nonhydro'
    
dirout = '/lustre/scratch4/turquoise/.mdt3/jrub/FIGS_HOMME/' 
moddir = '/lustre/scratch4/turquoise/.mdt3/jrub/Ertel/build/dcmip_tests/bw_l_'+dyn+'_'+str(X)+'radius_'+str(nlayer)+'level_ne'+str(ne)+'/theta/movies/'

ifile = 'dcmip2016_test11.nc'
outf = 'bw_'+dyn+'_'+str(X)+'X_'+str(nlayer)+'level_ne'+str(ne)
outf = outf +'_Ertel_lon'+str(lon0)


tframe = 0

ds = xr.open_dataset(moddir+ifile)
time = ds.time 
g = ds.geo.sel(time=time[tframe],lon=lon0) / (9.81 * 1000)


latmod = ds.lat 
z = np.tile(latmod,[g.shape[0],1])
f = 2 * 2 * np.pi * np.sin(np.pi * z / 180) * X / 86400


zeta_i = ds.zeta_x.sel(time=time[tframe],lon=lon0)
zeta_j = ds.zeta_y.sel(time=time[tframe],lon=lon0)
zeta_k = ds.zeta.sel(time=time[tframe],lon=lon0) #+ f

gradTh_i = ds.gradTh_x.sel(time=time[tframe],lon=lon0)
gradTh_j = ds.gradTh_y.sel(time=time[tframe],lon=lon0)
gradTh_k = ds.gradTh_z.sel(time=time[tframe],lon=lon0)

Va_gTh_i = zeta_i * gradTh_i
Va_gTh_j = zeta_j * gradTh_j
Va_gTh_k = zeta_k * gradTh_k
Va_gTh = Va_gTh_i + Va_gTh_j + Va_gTh_k

u   = ds.u.sel(time=time[tframe],lon=lon0)
Th  = ds.Th.sel(time=time[tframe],lon=lon0)
w   = ds.w.sel(time=time[tframe],lon=lon0)
rho = ds.rho.sel(time=time[tframe],lon=lon0)

#++++++++++++++++++
nk = 12
ml=0.1
mr=0.1
mt=0.05
mb=0.1
nx=3
ny=4

sx=0.05
sy=0.05
Fx=(1-(sx*(nx-1)+ml+mr))/nx
Fy=(1-(sy*(ny-1)+mt+mb))/ny


fig = plt.figure(figsize=(18,10.5), frameon=False) 

for iframe in np.arange(0,nk):
    R=int(np.ceil((iframe+1)/nx))
    C=int((iframe+1)-((R-1)*nx))
    Cx=(ml+(C-1)*Fx+(C-1)*sx)
    Cy=(1-mt-(Fy*R)-sy*(R-1))

    pos=(Cx,Cy,Fx,Fy)
    sfig=fig.add_axes(pos,frameon=True)


    v = eval(variables[iframe])

    cmap = eval('plt.cm.'+colmaps[iframe])
    cs = plt.contourf(z,g,v,np.linspace(limits[iframe][0],limits[iframe][1],50,endpoint=True),
                      cmap=cmap,linewidths=1.0) 
    #cs = plt.contourf(z,g,v,50,cmap=cmap,linewidths=1.0) 

    plt.colorbar()

    if (variables[iframe] == 'Th'):
        cs = plt.contour(z,g,v,np.linspace(limits[iframe][0],limits[iframe][1],30,endpoint=True),
                         linewidths=0.5,colors='black')

        #M =   f * ( (6357e3 / X ) * np.pi * z / 180) - u
        #sys.exit()
        #cs = plt.contour(z,g,M,30,linewidths=1.5,colors='white')

    plt.text(-89,35.2,label[iframe],fontsize=16,bbox=dict(facecolor='white',edgecolor='black',pad=2))


sys.exit()

#plt.savefig(figname)

