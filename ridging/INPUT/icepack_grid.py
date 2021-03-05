from midas.rectgrid import *
from midas.rectgrid_gen import *
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


xlen=5.;ys=60.;nxp=51
D0=500. # basin depth
xmid=0.5
ymid=0.5


if np.mod(nxp,20)==0: print('number of cell bounds must be uneen')
ncells=int((nxp-1)/2)


x=np.linspace(0.,1,nxp)
X,Y=np.meshgrid(x,x)
X=X*xlen
Y=Y*xlen+ys
sgrid=supergrid(xdat=X,ydat=Y,axis_units='degrees')
sgrid.grid_metrics()
grid=quadmesh(supergrid=sgrid)

#Construct the tracer grid
xh=X[1::2,1::2]
yh=xh.T
dpth=np.zeros((ncells,ncells))+D0
#coastline
xc=xh[0,:];yc=yh[:,0]
nc2=np.int(ncells / 2)
for k in np.arange(nc2):
    dpth[k,k+nc2:]=0.
    dpth[ncells-k-1,k+nc2:]=0.

#Save hgrid file
sgrid.write_nc('ocean_hgrid.nc')
#Save topog file
f=nc.Dataset('topog.nc','w')
f.createDimension('nx',ncells)
f.createDimension('ny',ncells)
f.createDimension('string',255)
tile=f.createVariable('tile','S1',('string'))
top=f.createVariable('depth','f8',('ny','nx'))
tile[:4]='tile1'
v=f.variables['tile']
d=np.empty(1,'S'+repr(len(v)))
d[0]='tile1'
dc=nc.stringtochar(d)
v[:]=dc
top[:]=dpth
f.close()


f=nc.Dataset('uvwind.nc','w')
f.createDimension('lon',ncells)
f.createDimension('lat',ncells)
f.createDimension('time',None)

lonv=f.createVariable('lon','f8',('lon',))
latv=f.createVariable('lat','f8',('lat',))
tv=f.createVariable('time','f8',('time',))
tv.units='days since 0001-01-01 00:00:00'
tv.calendar='julian'
#tv.modulo=' '
uv=f.createVariable('uwind','f8',('time','lat','lon'))
vv=f.createVariable('vwind','f8',('time','lat','lon'))
lonv[:]=xc
latv[:]=yc
V0=20.0;V0_len=1.0
U0=10.0;U_ref=1.0
v_wind = np.zeros(xh.shape)
time=0.
for i in np.arange(120):
    uv[i,:]=np.sin(i*np.pi/4.)*U0 + v_wind + U_ref
    vv[i,:]=np.sin(i*np.pi/4.)*V0 + v_wind
    tv[i]=time
    time=time+V0_len/8.
f.close()


f=nc.Dataset('tair.nc','w')
f.createDimension('lon',ncells)
f.createDimension('lat',ncells)
f.createDimension('time',None)

lonv=f.createVariable('lon','f8',('lon',))
latv=f.createVariable('lat','f8',('lat',))
tv=f.createVariable('time','f8',('time',))
tv.units='days since 0001-01-01 00:00:00'
tv.calendar='julian'
tv.modulo=' '
vv=f.createVariable('tair','f8',('time','lat','lon'))
lonv[:]=xc
latv[:]=yc
T_REF=273.15;T_RANGE=20.
tair = np.zeros(xh.shape)+T_REF
time=0.
for i in np.arange(120):
    vv[i,:]=np.sin(i*np.pi/4.)*T_RANGE + tair
    tv[i]=time
    time=time+V0_len/8.
f.close()

f=nc.Dataset('ice_ic.nc','w')
f.createDimension('lon',ncells)
f.createDimension('lat',ncells)
f.createDimension('time',None)

lonv=f.createVariable('lon','f8',('lon',))
latv=f.createVariable('lat','f8',('lat',))
tv=f.createVariable('time','f8',('time',))
tv.units='days since 0001-01-01 00:00:00'
tv.calendar='julian'
tv.modulo=' '
vv=f.createVariable('sit','f8',('time','lat','lon'))
cv=f.createVariable('sic','f8',('time','lat','lon'))
lonv[:]=xc
latv[:]=yc
sit = np.cos(4*(yh-ys)/xlen*np.pi)*2.0 + 3.0
sit[dpth==0]=0
sic = 0.5*np.cos(4*(yh-ys)/xlen*np.pi)+0.5
sic[dpth==0]=0
time=0.
for i in np.arange(1):
    vv[i,:]=sit
    cv[i,:]=sic
    tv[i]=time
    time=time+V0_len/8.
f.close()


######################
#Write Mosaics
######################

f=nc.Dataset('ocean_mosaic.nc','w')
f.createDimension('ntiles',1)
f.createDimension('string',255)
mc=f.createVariable('mosaic','c',('string',))
mc.standard_name = 'grid_mosaic_spec'
mc.children = 'contacts'
mc.grid_descriptor = ''
gridLoc = f.createVariable('gridlocation','c',('string',))
gridLoc.standard_name = 'grid_file_location'
gridFiles = f.createVariable('gridfiles','c',('ntiles','string',))
gridTiles = f.createVariable('gridtiles','c',('ntiles','string',))
f.grid_version = '0.2'
mc[:] = '\000' * 255
mc[:12] = 'ocean_mosaic'
gridLoc[:] = '\000' * 255
gridLoc[:2] = './'
gridFiles[:] = '\000' * 255
gridFiles[0,:14] = 'ocean_hgrid.nc'
gridTiles[:] = '\000' * 255
gridTiles[0,:5] = 'tile1'
f.close()

Ah=grid.Ah #Tracer cell area ('m2')
nl=len(numpy.where(dpth==0.)[0])
rg = nc.Dataset('atmos_mosaic_tile1Xland_mosaic_tile1.nc','w',format='NETCDF3_CLASSIC') # atmos_mosaic_tile1Xland_mosaic_tile1.nc
rg2 = nc.Dataset('land_mask.nc','w',format='NETCDF3_CLASSIC') # atmos_mosaic_tile1Xland_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('ncells',nl)  #It is unclear whether this works when nl=0. It does work for nl>0
rg.createDimension('two',2)
contact = rg.createVariable('contact','c',('string',))
contact.standard_name = 'grid_contact_spec'
contact.contact_type = 'exchange'
contact.parent1_cell = 'tile1_cell'
contact.parent2_cell = 'tile2_cell'
contact.xgrid_area_field = 'xgrid_area'
contact.distant_to_parent1_centroid = 'tile1_distance'
contact.distant_to_parent2_centroid = 'tile2_distance'
tile1_cell = rg.createVariable('tile1_cell','i4',('ncells','two',))
tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
tile2_cell = rg.createVariable('tile2_cell','i4',('ncells','two',))
tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
xgrid_area.standard_name = 'exchange_grid_area'
xgrid_area.units = 'm2'
tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
rg.grid_version = '0.2'

ni=ncells;nj=ncells
rg2.createDimension('nx',ni)
rg2.createDimension('ny',nj)
mask=rg2.createVariable('mask','f8',('ny','nx'))
mask.standard_name  = 'land fraction at T-cell centers'
mask.units = 'none'
mask[:,:]=0.0
rg2.grid_version = '0.2'

contact[:] = '\000' * 255
contact[:37] = 'atmos_mosaic:tile1::land_mosaic:tile1'
count=0
for j in range(nj):
    for i in range(ni):
        if dpth[j,i]==0.:
            tile1_cell[count] = [i+1,j+1]
            tile2_cell[count] = [i+1,j+1]
            tile1_distance[count] = [0.,0.]
            tile2_distance[count] = [0.,0.]
            xgrid_area[count] = Ah[j,i]
            count=count+1
            mask[j,i]=1.0
rg.close()
rg2.close()

rg = nc.Dataset('atmos_mosaic_tile1Xocean_mosaic_tile1.nc','w',format='NETCDF3_CLASSIC') # atmos_mosaic_tile1Xocean_mosaic_tile1.nc
rg2 = nc.Dataset('ocean_mask.nc','w',format='NETCDF3_CLASSIC') # atmos_mosaic_tile1Xland_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('ncells',ni*nj-nl) # -1 is for a single land point
print('ncells= ',ni*nj-nl)
rg.createDimension('two',2)
contact = rg.createVariable('contact','c',('string',))
contact.standard_name = 'grid_contact_spec'
contact.contact_type = 'exchange'
contact.parent1_cell = 'tile1_cell'
contact.parent2_cell = 'tile2_cell'
contact.xgrid_area_field = 'xgrid_area'
contact.distant_to_parent1_centroid = 'tile1_distance'
contact.distant_to_parent2_centroid = 'tile2_distance'
tile1_cell = rg.createVariable('tile1_cell','i4',('ncells','two',))
tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
tile2_cell = rg.createVariable('tile2_cell','i4',('ncells','two',))
tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
xgrid_area.standard_name = 'exchange_grid_area'
xgrid_area.units = 'm2'
tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
rg.grid_version = '0.2'
# Fill in data
contact[:] = '\000' * 255
contact[:38] = 'atmos_mosaic:tile1::ocean_mosaic:tile1'

rg2.createDimension('nx',ni)
rg2.createDimension('ny',nj)
mask=rg2.createVariable('mask','f8',('ny','nx'))
mask.standard_name  = 'ocean fraction at T-cell centers'
mask.units = 'none'
mask[:,:]=0.0
rg2.grid_version = '0.2'

count=0
for j in range(nj):
    for i in range(ni):
        if dpth[j,i]!=0:
            tile1_cell[count] = [i+1,j+1]
            tile2_cell[count] = [i+1,j+1]
            tile1_distance[count] = [0.,0.]
            tile2_distance[count] = [0.,0.]
            xgrid_area[count] = Ah[j,i]
            count=count+1
            mask[j,i]=1.0
rg.close()
rg2.close()

rg = nc.Dataset('land_mosaic_tile1Xocean_mosaic_tile1.nc','w',format='NETCDF3_CLASSIC') # land_mosaic_tile1Xocean_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('ncells',ni*nj-nl) # -1 is for a single land point
rg.createDimension('two',2)
contact = rg.createVariable('contact','c',('string',))
contact.standard_name = 'grid_contact_spec'
contact.contact_type = 'exchange'
contact.parent1_cell = 'tile1_cell'
contact.parent2_cell = 'tile2_cell'
contact.xgrid_area_field = 'xgrid_area'
contact.distant_to_parent1_centroid = 'tile1_distance'
contact.distant_to_parent2_centroid = 'tile2_distance'
tile1_cell = rg.createVariable('tile1_cell','i4',('ncells','two',))
tile1_cell.standard_name = 'parent_cell_indices_in_mosaic1'
tile2_cell = rg.createVariable('tile2_cell','i4',('ncells','two',))
tile2_cell.standard_name = 'parent_cell_indices_in_mosaic2'
xgrid_area = rg.createVariable('xgrid_area','f8',('ncells',))
xgrid_area.standard_name = 'exchange_grid_area'
xgrid_area.units = 'm2'
tile1_distance = rg.createVariable('tile1_distance','f8',('ncells','two'))
tile1_distance.standard_name = 'distance_from_parent1_cell_centroid'
tile2_distance = rg.createVariable('tile2_distance','f8',('ncells','two'))
tile2_distance.standard_name = 'distance_from_parent2_cell_centroid'
rg.grid_version = '0.2'
# Fill in data
contact[:] = '\000' * 255
contact[:37] = 'land_mosaic:tile1::ocean_mosaic:tile1'
count=0
for j in range(nj):
    for i in range(ni):
        if dpth[j,i]>0.:
            tile1_cell[count] = [i+1,j+1]
            tile2_cell[count] = [i+1,j+1]
            tile1_distance[count] = [0.,0.]
            tile2_distance[count] = [0.,0.]
            xgrid_area[count] = Ah[j,i]
            count=count+1


rg.close()

g = nc.Dataset('grid_spec.nc','w',format='NETCDF3_CLASSIC') # land_mosaic_tile1Xocean_mosaic_tile1.nc
rg.createDimension('string',255)
rg.createDimension('nfile_aXo',1) # -1 is for a single land point
rg.createDimension('nfile_aXl',1) # -1 is for a single land point
rg.createDimension('nfile_lXo',1) # -1 is for a single land point
atm_mosaic_dir = rg.createVariable('atm_mosaic_dir','c',('string',))
atm_mosaic_dir.standard_name = 'directory_storing_atmosphere_mosaic'
atm_mosaic_file = rg.createVariable('atm_mosaic_file','c',('string',))
atm_mosaic_file.standard_name = 'atmosphere_mosaic_file_name'
atm_mosaic = rg.createVariable('atm_mosaic','c',('string',))
atm_mosaic.standard_name = 'atmosphere_mosaic_name'
lnd_mosaic_dir = rg.createVariable('lnd_mosaic_dir','c',('string',))
lnd_mosaic_dir.standard_name = 'directory_storing_land_mosaic'
lnd_mosaic_file = rg.createVariable('lnd_mosaic_file','c',('string',))
lnd_mosaic_file.standard_name = 'land_mosaic_file_name'
lnd_mosaic = rg.createVariable('lnd_mosaic','c',('string',))
lnd_mosaic.standard_name = 'land_mosaic_name'
ocn_mosaic_dir = rg.createVariable('ocn_mosaic_dir','c',('string',))
ocn_mosaic_dir.standard_name = 'directory_storing_ocean_mosaic'
ocn_mosaic_file = rg.createVariable('ocn_mosaic_file','c',('string',))
ocn_mosaic_file.standard_name = 'ocean_mosaic_file_name'
ocn_mosaic = rg.createVariable('ocn_mosaic','c',('string',))
ocn_mosaic.standard_name = 'ocean_mosaic_name'
ocn_topog_dir = rg.createVariable('ocn_topog_dir','c',('string',))
ocn_mosaic_dir.standard_name = 'directory_storing_ocean_topog'
ocn_topog_file = rg.createVariable('ocn_topog_file','c',('string',))
ocn_topog_file.standard_name = 'ocean_topog_file_name'
aXo_file = rg.createVariable('aXo_file','c',('nfile_aXo','string',))
aXo_file.standard_name = 'atmXocn_exchange_grid_file'
aXl_file = rg.createVariable('aXl_file','c',('nfile_aXl','string',))
aXl_file.standard_name = 'atmXlnd_exchange_grid_file'
lXo_file = rg.createVariable('lXo_file','c',('nfile_lXo','string',))
lXo_file.standard_name = 'lndXocn_exchange_grid_file'

#Global attributes
rg.grid_version = '0.2'
rg.code_version = "$Name:  $"
rg.history = " "

atm_mosaic_dir[:] = '\000' * 255
atm_mosaic_dir[:2] = './'
atm_mosaic_file[:] = '\000' * 255
atm_mosaic_file[:15] = 'ocean_mosaic.nc'
atm_mosaic[:] = '\000' * 255
atm_mosaic[:12] = 'atmos_mosaic'
lnd_mosaic_dir[:] = '\000' * 255
lnd_mosaic_dir[:2] = './'
lnd_mosaic_file[:] = '\000' * 255
lnd_mosaic_file[:15] = 'ocean_mosaic.nc'
lnd_mosaic[:] = '\000' * 255
lnd_mosaic[:11] = 'land_mosaic'
ocn_mosaic_dir[:] = '\000' * 255
ocn_mosaic_dir[:2] = './'
ocn_mosaic_file[:] = '\000' * 255
ocn_mosaic_file[:15] = 'ocean_mosaic.nc'
ocn_mosaic[:] = '\000' * 255
ocn_mosaic[:12] = 'ocean_mosaic'
ocn_topog_dir[:] = '\000' * 255
ocn_topog_dir[:2] = './'
ocn_topog_file[:] = '\000' * 255
ocn_topog_file[:8] = 'topog.nc'
aXo_file[:,:] = '\000' * 255
aXo_file[:,:40] = 'atmos_mosaic_tile1Xocean_mosaic_tile1.nc'
aXl_file[:,:] = '\000' * 255
aXl_file[:,:39] = 'atmos_mosaic_tile1Xland_mosaic_tile1.nc'
lXo_file[:,:] = '\000' * 255
lXo_file[:,:39] = 'land_mosaic_tile1Xocean_mosaic_tile1.nc'

rg.close()

#S=state(grid=grid)
#vd={}
#vd['X']=ncells;vd['Y']=;vd['Z']=None;vd['T']=None
#vd['_FillValue']=-1.e20;vd['missing_value']=-1.e20
#S.add_field_from_array(dpth,'depth',var_dict=vd)

#S.write_nc('topog.nc',fields=['depth'])
