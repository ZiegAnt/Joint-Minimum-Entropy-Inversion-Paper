import numpy             as np
import matplotlib.pyplot as plt
import pygimli           as pg
from    pygimli.viewer.mpl      import createColorBarOnly


def plotJointResultsComparison(mesh, r_est_list, v_est_list,
                               c_ert, c_srt, lim_ert, lim_srt, label_list, labels=False,
                               ert_marks=None, srt_marks=None, ert_label=None, srt_label=None, marks=False,
                               figsize=(10,3.5)):
                               
    '''
    mesh: Inversion mesh (should be same for ERT and SRT)
    r_est_list/v_est_list: lists containing resistivity and velocity results (same order!)
    label_list: list containing labels that each row shares (only if labels=True)
    c_ert/srt: colormaps
    lim_ert/srt: colorbar limits
    ert/srt_marks: marks on respective colorbars with ert/srt_labels (only if marks=True)
    figsize: trivial
    '''
    
    if len(r_est_list)==1:
        fig, ax, [cax1, cax2] = plotResultsComparison(mesh, r_est_list[0], v_est_list[0], c_ert, c_srt, lim_ert, lim_srt,
                                                      ert_marks=ert_marks, srt_marks=srt_marks, 
                                                      ert_label=ert_label, srt_label=srt_label, marks=marks,
                                                      figsize=figsize)
    else:
        

        rows = len(r_est_list)
        fig, ax = plt.subplots(rows, 2, figsize=figsize)
        fig.tight_layout(pad=0.5)

        ax[0,0].set_title('ERT Results', fontsize = 16)
        ax[0,1].set_title('SRT Results', fontsize = 16)

        for i in range(rows): 
            pg.show(mesh, v_est_list[i], ax=ax[i,1], 
                        cMap=c_srt, cMin=lim_srt[0], cMax=lim_srt[1], 
                        colorBar=False, logScale=False)

            pg.show(mesh, r_est_list[i], ax=ax[i,0], 
                        cMap=c_ert, cMin=lim_ert[0], cMax=lim_ert[1], 
                        colorBar=False, logScale=True)

            if labels:
                for axis in [ax[i,0],ax[i,1]]:
                    axis.text(-32.5, -17, label_list[i], fontsize=10,
                              bbox={'facecolor': 'white', 'alpha': 0.8, 'pad': 5})

        # Adjust axis labels
        for axis in ax[:,0]:
            axis.set_ylabel('Z (m)', fontsize=12)
        for axis in ax[-1,:]:
            axis.set_xlabel('X (m)', fontsize=12)
        for ax_ar in ax[:-1,:]:
            for axis in ax_ar:
                axis.set_xticks([])
        for axis in ax[:,1]:
            axis.set_yticks([])
        for ax_ar in ax:
            for axis in ax_ar:
                axis.set_xlim(-35,35)
                axis.set_ylim(-20,0)
        for axis in ax:
            for a in axis:
                a.tick_params(labelsize=11)

        # Add colorbar ERT
        cax1 = ax[-1,0].inset_axes([0.05, -0.7, 0.9, 0.15])
        createColorBarOnly(ax=cax1, cMin=lim_ert[0], cMax=lim_ert[1], logScale=True,cMap=c_ert,
                          label=pg.unit('res'), orientation='horizontal')
        if marks:
            for i, v in enumerate(ert_marks):
                cax1.plot([v]*2, [0,1], 'w')
                cax1.text(v, 1.5, ert_label[i], fontsize=10, horizontalalignment='center', verticalalignment='center')

        # Add colorbar SRT
        cax2 = ax[-1,1].inset_axes([0.05, -0.7, 0.9, 0.15])
        createColorBarOnly(ax=cax2, cMin=lim_srt[0], cMax=lim_srt[1], logScale=False,cMap=c_srt,
                          label=pg.unit('vel'), orientation='horizontal')
        if marks:
            for i, r in enumerate(srt_marks):
                cax2.plot([r]*2, [0,1], 'w')
                cax2.text(r, 1.5, srt_label[i], fontsize=10, horizontalalignment='center', verticalalignment='center')

        # some axis edits
        cax1.xaxis.label.set_size(12)
        cax2.xaxis.label.set_size(12)
        cax1.tick_params(labelsize=11)
        cax2.tick_params(labelsize=11)

    return fig, ax, [cax1, cax2]


def plotResultsComparison(mesh, r_est, v_est,c_ert, c_srt, lim_ert, lim_srt,
                          ert_marks=None, srt_marks=None, ert_label=None, srt_label=None, marks=False,
                          figsize=(10,3.5)):
                               
    '''
    mesh: Inversion mesh (should be same for ERT and SRT)
    r_est/v_est: resistivity and velocity results
    c_ert/srt: colormaps
    lim_ert/srt: colorbar limits
    ert/srt_marks: marks on respective colorbars with ert/srt_labels (only if marks=True)
    figsize: trivial
    '''

    
    fig, ax = plt.subplots(1, 2, figsize=figsize)
    fig.tight_layout(pad=0.5)

    ax[0].set_title('ERT Results', fontsize = 14)
    ax[1].set_title('SRT Results', fontsize = 14)

    pg.show(mesh, v_est, ax=ax[1], 
                    cMap=c_srt, cMin=lim_srt[0], cMax=lim_srt[1], 
                    colorBar=False, logScale=False)

    pg.show(mesh, r_est, ax=ax[0], 
                    cMap=c_ert, cMin=lim_ert[0], cMax=lim_ert[1], 
                    colorBar=False, logScale=True)


    # Adjust axis labels
    ax[0].set_ylabel('Z (m)')
    ax[0].set_xlabel('X (m)')
    ax[1].set_xlabel('X (m)')
    ax[1].set_yticks([])
    for axis in ax:
            axis.set_xlim(-35,35)
            axis.set_ylim(-20,0)

    # Add colorbar ERT
    cax1 = ax[0].inset_axes([0.05, -0.6, 0.9, 0.15])
    createColorBarOnly(ax=cax1, cMin=lim_ert[0], cMax=lim_ert[1], logScale=True,cMap=c_ert,
                      label=pg.unit('res'), orientation='horizontal')
    if marks:
        for i, v in enumerate(ert_marks):
            cax1.plot([v]*2, [0,1], 'w')
            cax1.text(v, 1.3, ert_label[i], fontsize=8, horizontalalignment='center', verticalalignment='center')

    # Add colorbar SRT
    cax2 = ax[1].inset_axes([0.05, -0.6, 0.9, 0.15])
    createColorBarOnly(ax=cax2, cMin=lim_srt[0], cMax=lim_srt[1], logScale=False,cMap=c_srt,
                      label=pg.unit('vel'), orientation='horizontal')
    if marks:
        for i, r in enumerate(srt_marks):
            cax2.plot([r]*2, [0,1], 'w')
            cax2.text(r, 1.3, srt_label[i], fontsize=8, horizontalalignment='center', verticalalignment='center')
        
    return fig, ax, [cax1, cax2]

def showPseudosections(geoContainer, rhoa, idx_list, Lines, ax, clim, cmap='gnuplot', annotation=True, Type='data', colorBar=True):
    '''
    geoContainer: dataContainer containing a,b,m,n
    rhoa: numpy array containing apparent resistivities or misfits
    idx_list: List containing [0, end(line1), ...]
    Lines: number of lines that are present in the data
    ax: matplotlib subplot axis to use for plotting
    clim: [cmin, cmax]
    cmap: colormap
    annotation: if True then Line description is added to pseudosections
    Type: either 'data' or 'misfit'
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
                             
    ax.set_xlabel('Midpoint index')
    ax.set_ylabel('Level')
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
    
    resp = mgr.inv.response.array()  #get response
    
    # apparent velocities
    va_obs = offset/data
    va_pre = offset/resp
    
    # Calc apparent velocity misfit in %
    misfit =  100*(va_obs-va_pre)/(va_obs)
    
    #Transfer Vectors into 2D Matrix
    misfit_mat = np.zeros((mgr.data.sensorCount(),mgr.data.sensorCount()))
    for i in range(mgr.data.sensorCount()):
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
    ylab = [item.get_text() for item in ax.get_yticklabels()]
    ylab[1:] = [f'{d_shot*int(y)}' for y in ylab[1:]]
    ax.set_yticklabels(ylab)
    
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