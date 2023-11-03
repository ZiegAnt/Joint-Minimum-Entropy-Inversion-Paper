import numpy             as np
import matplotlib.pyplot as plt
import pygimli           as pg
from pygimli.viewer.pv import pgMesh2pvMesh  
from    pygimli.viewer.mpl      import createColorBarOnly
import pyvista as pyv
from pygimli.core import getRotation, degToRad
import pygimli.meshtools as mt
from pygimli.physics.gravimetry import MagneticsModelling

def showPseudosections(geoContainer, rhoa, idx_list, Lines, ax, clim, cmap='gnuplot', annotation=True, Type='data', colorBar=True):
    '''
    geoContainer: dataContainer containing a,b,m,n
    rhoa: numpy array containing apparent resistivities
    Lines: number of lines that are present in the data
    ax: matplotlib subplot axis to use for plotting
    clim: [cmin, cmax]
    cmap: colormap
    annotation: if True then Line description is added to pseudosections
    Type: either 'data' or 'misfit'
    colorBar: if True than with colorbar if False then  not
    '''
    x = []; y = []; c = []; s = []
    
    a = geoContainer['a']
    b = geoContainer['b']
    m = geoContainer['m']
    n = geoContainer['n']

    mp = 0.5*(a+m) #midpoint
    diplen = b-a # dipole length
    lev = (m-b) # dipole separation/level

    mdl = max(diplen) # maximum dipole length 
    dl_shift = 1/(mdl/2) # small shift for dipole length
    level_new = lev + dl_shift*diplen # resulting y-value
        
    d = rhoa
    
    x_new = mp.tolist()
    y_new = level_new.tolist()
    c_new = d.tolist()
    s_new = diplen.tolist()
    
    if Type=='data':
        kw = {'cmap': cmap, 'norm': 'log', 'vmin':clim[0], 'vmax':clim[1]}
        
        for i in range(Lines):
            im = ax.tricontourf(x_new[idx_list[i]:idx_list[i+1]],
                                y_new[idx_list[i]:idx_list[i+1]],
                                c_new[idx_list[i]:idx_list[i+1]], 
                                levels=np.arange(clim[0], clim[1], 1), extend='both',**kw)
        im.cmap.set_under('k')
        im.cmap.set_over('y')
    
    if Type=='misfit':
        kw = {'cmap': cmap, 'vmin':clim[0], 'vmax':clim[1]}
        for i in range(Lines):
            im = ax.tricontourf(x_new[idx_list[i]:idx_list[i+1]],
                                y_new[idx_list[i]:idx_list[i+1]],
                                c_new[idx_list[i]:idx_list[i+1]], 
                                levels=np.arange(clim[0], clim[1], 1), extend='both',**kw)
        im.cmap.set_under('b')
        im.cmap.set_over('r')
                             
    ax.set_xlabel('Midpoint index',fontsize=10)
    ax.set_ylabel('Level',fontsize=10)
    ax.invert_yaxis()
    
    if annotation:
        texty = max(y_new)-5
        for i in range(Lines):
            textx = i*max(x_new)/(Lines)+15
            ax.text(textx, texty, f'Line {i+1}', horizontalalignment='center', verticalalignment='center')
            
    if colorBar:
        cax = ax.inset_axes([1.02, 0.05, 0.015, 0.9])    
        if Type=='data':
            createColorBarOnly(ax=cax, cMin=clim[0], cMax=clim[1], logScale=True,cMap=cmap,
                              label=pg.unit('rhoa'), orientation='vertical')
        else:
            createColorBarOnly(ax=cax, cMin=clim[0], cMax=clim[1], logScale=True,cMap=cmap,
                              label='Misfit (%)', orientation='vertical')
        return ax, cax
    else:
        return ax

def showTTMatrix(mgr, Type, ax, lim=None, cmap='jet', colorBar=True): 
    '''
    mgr: traveltime manager
    Type: either 'misfit', 'appvel_pre' or 'appvel_obs'
    ax: plt.mappable
    lim: limits [min,max]
    cmap: colormap
    colorBar: Plot colorbar right to Figure or not
    '''
    # get data from manager
    shot = mgr.data['s']
    rec = mgr.data['g']
    data = mgr.data['t'].array()
    offset = [np.linalg.norm(mgr.data.sensorPosition(shot[i]).array()-mgr.data.sensorPosition(rec[i]).array()) for i in range(len(shot))]
    d_shot = np.unique(shot)[1]-np.unique(shot)[0] # shot index separation
    
    resp = mgr.inv.response.array()  #et response
    
    # apparent velocities
    va_obs = offset/data
    va_pre = offset/resp
    
    # Calc apparent velocity misfit in %
    misfit =  100*(va_obs-va_pre)/(va_obs)
    
    #Transfer Vectors into 2D Matrix
    misfit_mat = np.zeros((data.sensorCount(),data.sensorCount()))
    for i in range(data.sensorCount()):
        misfit_mat[i,i] = None
    app_vel_obs = misfit_mat.copy()
    app_vel_pre = misfit_mat.copy()
    for i in range(len(shot)):
        misfit_mat[shot[i],rec[i]] = misfit[i]
        app_vel_obs[shot[i],rec[i]] = va_obs[i]
        app_vel_pre[shot[i],rec[i]] = va_pre[i]
    misfit_mat = misfit_mat[::2,:]
    app_vel_obs = app_vel_obs[::2,:]
    app_vel_pre = app_vel_pre[::2,:]

    # Plot 2D Matrix
    if lim==None:
        if Type=='misfit':
            im = ax.imshow(misfit_mat, cmap=cmap)
            lim = [min(misfit),max(misfit)]
        if Type=='appvel_obs':
            im = ax.imshow(app_vel_obs, cmap=cmap)
            lim = [min(va_obs),max(va_obs)]
        if Type=='appvel_pre':
            im = ax.imshow(app_vel_pre, cmap=cmap)
            lim = [min(va_pre),max(va_pre)]
    else:
        if Type=='misfit':
            im = ax.imshow(misfit_mat, cmap=cmap, vmin=lim[0], vmax=lim[1])
        if Type=='appvel_obs':
            im = ax.imshow(app_vel_obs, cmap=cmap, vmin=lim[0], vmax=lim[1])
        if Type=='appvel_pre':
            im = ax.imshow(app_vel_pre, cmap=cmap, vmin=lim[0], vmax=lim[1])
    # Axes labels
    ax.set_xlabel('Sensor index')
    ax.set_ylabel('Source index')
    
    # Colorbar
    if colorBar:
        cax = ax.inset_axes([1.02, 0.05, 0.015, 0.9])
        if Type=='misfit':
            createColorBarOnly(ax=cax, cMin=lim[0], cMax=lim[1], logScale=False,cMap=cmap,
                               label='Misfit (%)', orientation='vertical')
        else:
            createColorBarOnly(ax=cax, cMin=lim[0], cMax=lim[1], logScale=False,cMap=cmap,
                               label=pg.unit('va'), orientation='vertical')
        return ax, cax
    else:
        return ax


def showMagResult(pnts, inv):
    d_obs = inv.dataVals
    d_pre = inv.response
    misfit = 100*(d_obs-d_pre)/d_obs          # Misfit in %
    
    fig, ax = plt.subplots(ncols=3, figsize=(15, 4.5), sharex=True, sharey=True, constrained_layout=True)
    mm = np.max(np.abs(d_obs))
    
    im0 = ax[0].scatter(pnts[:,0], pnts[:,1], c=d_obs, marker='8', cmap="PuOr", vmin=-mm, vmax=mm)
    im1 = ax[1].scatter(pnts[:,0], pnts[:,1], c=d_pre, marker='8', cmap="PuOr", vmin=-mm, vmax=mm)
    im2 = ax[2].scatter(pnts[:,0], pnts[:,1], c=misfit, marker='8', cmap="coolwarm", vmin=-10, vmax=10)
    
    ax[0].set_title('Measured Data')
    ax[1].set_title('Model Response')
    ax[2].set_title(f'Misfit $\chi^2$ = {inv.chi2():.2f}')
    ax[0].set_ylabel('Northing in m')
    ax[0].set_xlabel('Easting in m')
    ax[1].set_xlabel('Easting in m')
    ax[2].set_xlabel('Easting in m')
    cb1 = fig.colorbar(im0, ax=ax[0], orientation='vertical')
    cb2 = fig.colorbar(im1, ax=ax[1], orientation='vertical')
    cb3 = fig.colorbar(im2, ax=ax[2], orientation='vertical')

    cb1.ax.set_title('TFA\nnT')
    cb2.ax.set_title('TFA\nnT')
    cb3.ax.set_title('RMS\n%')
    return fig, ax
    
def showMagMisfit(pnts, input_type, data_input, ax, cmap='coolwarm', lim=20, size=3):
    if input_type=='Inversion_framework':
         inv = data_input
         d_obs = inv.dataVals
         d_pre = inv.response
         misfit = 100*(d_obs-d_pre)/d_obs          # Misfit in %
    else:
        misfit = data_input
        

    ax.scatter(pnts[:,0], pnts[:,1], c=misfit, s=size, marker='8', cmap="coolwarm", vmin=-lim, vmax=lim)
    ax.set_ylabel('Northing in m',fontsize=10)
    ax.set_xlabel('Easting in m',fontsize=10)

    return ax
     
  
def AZdrawSlice(ax, mesh,cmap, normal=[1, 0, 0], **kwargs):
    label = kwargs.pop('label', None)
    data = kwargs.pop('data', None)
    mesh = pgMesh2pvMesh(mesh, data, label)

    try:
        single_slice = mesh.slice(normal, **kwargs)

    except AssertionError as e:
        # 'contour' kwarg only works with point data and breaks execution
        pg.error(e)
    else:
        # REVIEW: bounds and axes might be confused with the outline..?!
        outline = mesh.outline()
        ax.add_mesh(outline, color="k")
        ax.add_mesh(single_slice, cmap=cmap)

    return ax

def AZdrawSlice_along_line(ax, mesh, cmap, pnts, **kwargs):
    '''
    example: AZdrawSlice_along_line(pl, mesh, cmap=c_mag, pnts=p1, data=mesh["sus"], label="sus")

    Parameters
    ----------
    ax : pyvista.plotter
    mesh : 3D Mesh
    cmap : colormap (e.g. 'jet')
    pnts : np.array() containing 3D coordinates for points along the line
    **kwargs : for add_mesh command
    '''    
    label = kwargs.pop('label', None)
    data = kwargs.pop('data', None)
    mesh = pgMesh2pvMesh(mesh, data, label)

    sec = pyv.Spline(pnts)
    
    slc = mesh.slice_along_line(sec)
    ax.add_mesh(slc, cmap=cmap)
    
    # outline = mesh.outline()
    # pl.add_mesh(mesh.outline(), color="k")
    

    return ax

def interpolate3D_to_2DSection(mesh_list, point_list, mesh3D, data3D):
    '''
    interpolate 3D data to 2D section with linear interpolation
    '''

    data_int = []

    for i,p in enumerate(point_list):
        m = mesh_list[i]


        start = pg.RVector3((p[0][0],p[0][1],0))
        end = pg.RVector3((p[-1][0],p[-1][1],0))

        src = pg.RVector3(0.0, 0.0, 0.0).norm(pg.RVector3(0.0, 0.0, -10.0),pg.RVector3(10.0, 0.0, -10.0))
        dest = start.norm(start - pg.RVector3(0.0, 0.0, 10.0), end)
        rot = getRotation(src, dest)


        m.rotate(degToRad(pg.RVector3(90.0, 0.0, 0.0))) # puts hgeight from y to z coordinate
        m.transform(rot)                                # rotates it correctly
        m.translate(start)                              # moves (x,y)=(0,0) [i.e. line start] to correct point

        d_int = mt.interpolate(m, mesh3D, data3D)
        data_int.append(d_int)

        # Reverse translations and rotations to initial locations
        m.translate(-start)
        m.transform(pg.core.RMatrix(np.array([1,-1,1])*np.array(rot)))
        m.rotate(degToRad(pg.RVector3(-90.0, 0.0, 0.0)))

    return data_int

def Move2DSectionTo3D(mesh_list, point_list):
    '''
    interpolate 3D data to 2D section with linear interpolation
    '''
    mesh_list3D = []

    for i,p in enumerate(point_list):
        m = mesh_list[i]


        start = pg.RVector3((p[0][0],p[0][1],0))
        end = pg.RVector3((p[-1][0],p[-1][1],0))

        src = pg.RVector3(0.0, 0.0, 0.0).norm(pg.RVector3(0.0, 0.0, -10.0),pg.RVector3(10.0, 0.0, -10.0))
        dest = start.norm(start - pg.RVector3(0.0, 0.0, 10.0), end)
        rot = getRotation(src, dest)
        
        m.rotate(degToRad(pg.RVector3(90.0, 0.0, 0.0))) # puts hgeight from y to z coordinate
        m.transform(rot)                                # rotates it correctly
        m.translate(start)                              # moves (x,y)=(0,0) [i.e. line start] to correct point

        mesh_list3D.append(m)
    return mesh_list3D


def createMagManager(mesh, df):
    '''
    mesh: Inversion mesh
    df: pandas dataframe holding all magnetic data (including geometry)
    
    Returns: Magnetics manager and magnetic data in pg.Vector format
    '''
    # PARAMETERS (based on Boxberg 2011)
    F = 48487.4   # Median field intensity in nT
    I = 65.70     # Median inclination
    D = 0.85      # Median declination
    
    # Get data
    d_obs = df['F'].to_numpy()-F   # Total Field Anomaly numpy array
    d_mag = pg.Vector(d_obs)       # Total Field Anomaly RVector
    noise = df['noise'].to_numpy() # Estimated noise of measurements (by device)
    pnts = np.array([[df['X'][i], df['Y'][i], df['Z'][i]] for i in range(len(df))]) # Sensor locations
    
    # Natural magnetic field
    H = F * np.cos(np.deg2rad(I))
    X = H * np.cos(np.deg2rad(D))
    Y = H * np.sin(np.deg2rad(D))
    Z = F * np.sin(np.deg2rad(I))
    igrf = [D, I, H, X, Y, Z, F]

    # Foraward Operator
    cmp = ["TFA"]  # Total Field Anomaly
    fop = MagneticsModelling(mesh, pnts, cmp, igrf)
    
    # Set relError
    relError = np.zeros(len(df)) # calculate relative error
    for i in range(len(df)):
        relError[i] = abs(noise[i]/d_mag[i])+0.02 # additional 2% rel error
        
    # Set Inversion FW
    inv = pg.Inversion(fop=fop, verbose=True)
    inv.setRegularization(limits=[0, 0.15])  # to limit values, diatreme should have around 0.1
    
    #Set manager
    MAG = pg.frameworks.methodManager.MethodManager(fop=fop, fw=inv, data=d_mag)
    
    return MAG, d_mag, relError
    
    
from pygimli.utils import ProgressBar

def interpolate2D_to_3DIndexNN(mesh_list, point_list, mesh3D):
    '''
    Get indices that can be used to fill 2D section parameters wrt 3D mesh using nearest Neighbour interpolation.
    '''
    int_idx_all = []
    pBar_mesh = ProgressBar(its=len(mesh_list), width=40, sign='+')

    for i,p in enumerate(point_list):
        m = mesh_list[i]
        int_idx = []
                            
        start = pg.RVector3((p[0][0],p[0][1],0))
        end = pg.RVector3((p[-1][0],p[-1][1],0))

        src = pg.RVector3(0.0, 0.0, 0.0).norm(pg.RVector3(0.0, 0.0, -10.0),pg.RVector3(10.0, 0.0, -10.0))
        dest = start.norm(start - pg.RVector3(0.0, 0.0, 10.0), end)
        rot = getRotation(src, dest)


        m.rotate(degToRad(pg.RVector3(90.0, 0.0, 0.0))) # puts hgeight from y to z coordinate
        m.transform(rot)                                # rotates it correctly
        m.translate(start)                              # moves (x,y)=(0,0) [i.e. line start] to correct point

        pBar = ProgressBar(its=m.cellCount(), width=40, sign='+')
        for j,pos in enumerate(m.cellCenters()):
            cell = mesh3D.findCell(pos)
            if cell:
                int_idx.append(cell.id())
            else:
                int_idx.append('nan')

            pBar.update(j)

        # Reverse translations and rotations to initial locations
        m.translate(-start)
        m.transform(pg.core.RMatrix(np.array([1,-1,1])*np.array(rot)))
        m.rotate(degToRad(pg.RVector3(-90.0, 0.0, 0.0)))
                            
        pBar_mesh.update(i)
        int_idx_all.append(int_idx)

    return int_idx_all

def plotJointResultsComparison(mesh, coverage, res_list, sus_list,
                               c_ert, c_mag, lim_ert, lim_mag, label_list=[], labels=False,figsize=(10,3.5)):
                               
    
    if len(res_list)==1:
        fig, ax, [cax1, cax2] = plotResultsComparison(mesh,coverage,res_list[0], sus_list[0], c_ert, c_mag, lim_ert, lim_mag, figsize=figsize)
    else:
        rows = len(res_list)
        fig, ax = plt.subplots(rows, 2, figsize=figsize)
        fig.tight_layout(pad=0.1)

        ax[0,0].set_title('ERT Results', fontsize = 14)
        ax[0,1].set_title('MAG Results', fontsize = 14)

        for i in range(rows): 
            pg.show(mesh, sus_list[i], ax=ax[i,1], coverage=coverage,
                        cMap=c_mag, cMin=lim_mag[0], cMax=lim_mag[1], 
                        colorBar=False, logScale=False)

            pg.show(mesh, res_list[i], ax=ax[i,0], coverage=coverage,
                        cMap=c_ert, cMin=lim_ert[0], cMax=lim_ert[1], 
                        colorBar=False, logScale=True)

            if labels:
                for axis in [ax[i,0],ax[i,1]]:
                    axis.text(1, 414, label_list[i], fontsize=10,
                              bbox={'facecolor': 'white', 'alpha': 0.8, 'pad': 5})

        # Adjust axis labels
        for axis in ax[:,0]:
            axis.set_ylabel('Elevation (m.a.s.l.)')
        for axis in ax[-1,:]:
            axis.set_xlabel('Distance (m)')
        for ax_ar in ax[:-1,:]:
            for axis in ax_ar:
                axis.set_xticks([])
        for axis in ax[:,1]:
            axis.set_yticks([])
        for ax_ar in ax:
            for axis in ax_ar:
                axis.set_ylim(408,455)

        # Add colorbar ERT
        cax1 = ax[-1,0].inset_axes([0.05, -0.7, 0.9, 0.15])
        createColorBarOnly(ax=cax1, cMin=lim_ert[0], cMax=lim_ert[1], logScale=True,cMap=c_ert,
                          label=pg.unit('res'), orientation='horizontal')

        # Add colorbar MAG
        cax2 = ax[-1,1].inset_axes([0.05, -0.7, 0.9, 0.15])
        createColorBarOnly(ax=cax2, cMin=lim_mag[0], cMax=lim_mag[1], logScale=False,cMap=c_mag,
                          label='Magnetic Susceptibility', orientation='horizontal')

    return fig, ax, [cax1, cax2]


def plotResultsComparison(mesh, coverage, r_est, v_est,c_ert, c_mag, lim_ert, lim_mag,figsize=(10,3.5)):

    
    fig, ax = plt.subplots(1, 2, figsize=figsize)
    fig.tight_layout(pad=0.1)

    ax[0].set_title('ERT Results', fontsize = 14)
    ax[1].set_title('MAG Results', fontsize = 14)

    pg.show(mesh, v_est, ax=ax[1], coverage=coverage,
                    cMap=c_mag, cMin=lim_mag[0], cMax=lim_mag[1], 
                    colorBar=False, logScale=False)

    pg.show(mesh, r_est, ax=ax[0], coverage=coverage,
                    cMap=c_ert, cMin=lim_ert[0], cMax=lim_ert[1], 
                    colorBar=False, logScale=True)


    # Adjust axis labels
    ax[0].set_ylabel('Elevation (m.a.s.l.)')
    ax[0].set_xlabel('Distance (m)')
    ax[1].set_xlabel('Distance (m)')
    ax[1].set_yticks([])
    for axis in ax:
            axis.set_ylim(408,455)

    # Add colorbar ERT
    cax1 = ax[0].inset_axes([0.05, -0.7, 0.9, 0.15])
    createColorBarOnly(ax=cax1, cMin=lim_ert[0], cMax=lim_ert[1], logScale=True,cMap=c_ert,
                      label=pg.unit('res'), orientation='horizontal')

    # Add colorbar MAG
    cax2 = ax[1].inset_axes([0.05, -0.7, 0.9, 0.15])
    createColorBarOnly(ax=cax2, cMin=lim_mag[0], cMax=lim_mag[1], logScale=False,cMap=c_mag,
                      label='Magnetic Susceptibility', orientation='horizontal')
    
        
    return fig, ax, [cax1, cax2]