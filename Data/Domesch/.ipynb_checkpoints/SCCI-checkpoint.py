import os
import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


class SCCI(object):
    def __init__(self, managers=[], **kwargs):
        """Structurally coupled cooperative inversion framework."""
        # scci parameters
        self.a = None
        self.b = None
        self.c = None
        self.cmin = None
        self.cmax = None
        self.setCouplingPara()

        # init managers list
        self.managers = managers
        # I think we don't need the managers unless we plot etc.

        self.resultdir = None  # "scci_results"

        self.names = kwargs.pop("names",
                                ["M{0}".format(i+1) for i in range(10)])
        # init constraint weights

    def setCouplingPara(self, a=0.1, b=0.1, c=1.0, cmin=0.2, cmax=1.0):
        """
        Set coupling parameter.
        """
        self.a = a
        self.b = b
        self.c = c
        self.cmin = cmin
        self.cmax = cmax

    @property
    def invs(self):
        return self._gather('inv', raise_error=True)

    @property
    def fops(self):
        return self._gather('fop', raise_error=True)

    @property
    def npara(self):
        npara = 0
        paras = self._gather('_numPara', default=1)
        npara += np.sum(paras)
        return npara

    @property
    def nparas(self):
        paras = self._gather('_numPara', default=1)
        return paras

    @property
    def roughnesses(self):
        output = []
        # do for all managers == inversion instances
        for im, man in enumerate(self.managers):
            model = man.inv.model

            if man.fop.constraints().cols() != model.size():
                man.fop.createConstraints()

            roughness = man.inv.inv.pureRoughness(model)

            seg = int(len(roughness)/self.nparas[im])
            # do for all parameters
            for i in range(self.nparas[im]):
                output.append(np.array(roughness[i * seg:i * seg + seg]))

        return output

    def singleCWeights(self):
        cweights = []
        for rough in self.roughnesses:
            cweight = roughness2CWeight(
                rough, a=self.a, b=self.b, c=self.c, Min=0,  # self.cmin,
                Max=self.cmax)
            pg.debug([np.min(cweight), np.max(cweight)])
            cweights.append(cweight)

        return cweights

    def updateConstraintWeights(self):
        """
        Set new constraint weights based on the actual model roughnesses.
        """
        pg.debug('SCCI: updateConstraintWeights()')
        # write down single c weights
        single_weights = np.array(self.singleCWeights())
        pg.debug('single c weights before updateConstraintWeights: {}'
                 .format(single_weights))

        new_cweight = np.ones_like(single_weights)

        # each new c weight consists of the combination of cweights of all
        # other parameters
        for ipar in range(self.npara):
            all_others = np.arange(self.npara) != ipar
            pg.debug('SCCI: all others: {}'.format(all_others))
            new_cweight[all_others] *= single_weights[ipar]
            # new_cweight *= single_weights[ipar]

        # cut extreme values according to cmin, cmax
        new_cweight = np.minimum(new_cweight, self.cmax)
        new_cweight = np.maximum(new_cweight, self.cmin)

        pg.debug('min/max new c weight: {}/{}'
                 .format(np.min(new_cweight), np.max(new_cweight)))

        # set c weights
        total = 0
        matrices = []
        for iinv, inv in enumerate(self.invs):
            weight = []
            for j in range(self.nparas[iinv]):
                weight.extend(new_cweight[total])
                total += 1

            pg.debug('SCCI: manager {}, weights {}'.format(
                iinv, np.shape(weight)))

            inv.inv.setCWeight(weight)

            matrices.append(inv.fop.constraints())

        return new_cweight, matrices

    def runCoupled(self, save=False, **kwargs):
        """Run coupled inversion."""
        maxIter = kwargs.pop("maxIter", 8)
        if self.resultdir is None:
            self.resultdir = "result_a{}_b{}".format(self.a, self.b).replace(
                ".", "_")

        if save and not os.path.exists(self.resultdir):
            os.makedirs(self.resultdir)

        for i in range(maxIter):
            print("Coupled inversion {0}".format(i+1))
            self.updateConstraintWeights()
            for i, inv in enumerate(self.invs):
                basename = self.resultdir + "/" + self.names[i]
                if save:
                    np.save(basename+'_cWeight_{}.npy'.format(i + 1),
                            inv.inv.cWeight())

                inv.inv.oneStep()
                if save:
                    np.save(basename + '_response_{}.npy'.format(i + 1),
                            inv.response)
                    np.save(basename + '_model_{}.npy'.format(i + 1),
                            inv.model)

    def _gather(self, attribute, raise_error=False, default=None):
        """ Internal function.
        Gather variables from underlaying managers for convenience. """
        ret = []
        for mi, manager in enumerate(self.managers):
            if hasattr(manager, attribute):
                ret.append(getattr(manager, attribute))
            else:
                if raise_error:
                    raise AttributeError(
                        'Manager {} of type {} has no attribute "{}"'
                        .format(mi, type(manager), attribute))
                else:
                    ret.append(default)
        return ret

    def _gather_from_cpp(self, call, raise_error=False, default=None):
        """ Internal function.
        Gather variables from given callable from underlaying managers
        for convenience.
        """
        ret = []
        for mi, manager in enumerate(self.managers):
            if hasattr(manager, call):
                to_call = getattr(manager, call)
                if not callable(to_call):
                    if raise_error:
                        raise TypeError('{} object of manager {} of type {} '
                                        'is not callable'.format(type(to_call),
                                                                 mi,
                                                                 type(manager)))
                    else:
                        ret.append(default)
                else:
                    ret.append(to_call())
            else:
                if raise_error:
                    raise AttributeError(
                        'Manager {} of type {} has no function "{}"'
                        .format(mi, type(manager), call))
                else:
                    ret.append(default)
        return ret


def roughness2CWeight(vec, a=0.1, b=0.1, c=1.0, Max=1.0, Min=0.2):
    """
    Structural constraint weight as function of roughness. See
    Günther et al. (2010, SAGEEP extended abstract) (case a>0) for details.
    """
    avec = np.absolute(vec)
    cfun = (a / (avec + a) + b)**c
    # confine between min and max
    # cfun = (cfun - min(cfun)) / (max(cfun) - min(cfun)) * (Max-Min) + Min
    cfun = np.minimum(cfun, Max)
    cfun = np.maximum(cfun, Min)
    return cfun


def getNormalizedRoughness(mesh, model, trans='log', tmin=None, tmax=None):
    """ Returns normalized roughness of the model with respect to the given
    mesh and transformation.

    The roughness is defined on each cell boundary with connection to another
    cell (so not for cells at the boundary of the mesh), so the number of
    values is not number of nodes nor number of cells or boundaries and
    strongly depends on the mesh.
    """
    # Step 1: init of mesh, fop and inv
    # we just need the methods of the inv and fop instances, so no input needed
    mesh.createNeighbourInfos()

    fop = pg.core.ModellingBase()
    fop.setMesh(mesh)

    # pseudo data, fop, verbose, debug
    inv = pg.core.Inversion([1, 2, 3], fop, True, False)

    # Step 2: take care of transformation
    if trans.lower() == 'log':
        if tmin is not None and tmax is not None:
            transmodel = pg.core.TransLogLU(tmin, tmax)
        else:
            transmodel = pg.core.TransLogLU()

    elif trans.lower() == 'cot':
        if tmin is None:
            tmin = 0

        if tmax is None:
            tmax = 0.7

        transmodel = pg.core.TransCotLU(tmin, tmax)

    else:
        raise Exception(
            'Transformation not implemented, expect "log" or "cot", not "{}".'
            .format(trans))

    # Step 3: region manager initialization (enables use of createConstraints)
    inv.setTransModel(transmodel)
    fop.regionManager()
    fop.createConstraints()

    # Step 4: with active constraints for the mesh, we simply get the roughness
    roughness = inv.pureRoughness(model)

    # return as numpy array
    return np.asarray(roughness)

import numpy as np
# import pygimli as pg
import matplotlib.pyplot as plt
# import matplotlib.patches as mpatches


def lineWidthFromC(cw, cmin=0.1, cmax=1, lmin=0, lmax=0.8):
    """Compute linewidth from given c value given cmin/cmax and lmin/lmax."""
    linew = (1 - (cw - cmin) / (cmax - cmin))*(lmax - lmin) + lmin
    linew = np.max([np.min([linew, lmax]), lmin])
    return linew


def drawCWeight(ax, mesh, cweight, lmin=0, lmax=0.8, cmin=0.1, cmax=1,
                min_plot=0.02, color='black', cell_indices=None):
    """ Draws the given cweights defined for given mesh on given ax.

    Parameters
    ----------

    ax: plt.ax
        Ax to plot constraint weights in.

    mesh: pg.Mesh
        Mesh object where the constraints are defined in.

    cweight: np.ndarray
        Constraint values to be plotted.

    lmin: float [ 0 ]
        Minimum linewidth for maximum cweight defined via **cmax**.
        Note that by default high constraint values are plotted with thinner
        lines.

    lmin: float [ 0.8 ]
        Maximum linewidth for minimum cweight defined via **cmin**.

    cmin: float [ 0 ]
        Minimum constraint weight to plot. All values smaller than cmin are
        plotted with the same linewidth as cmin.

    cmax: float [ 1 ]
        Maximum constraint weight to plot. All values greater than cmax are
        plotted with the same linewidth as cmax.

    min_plot: float [ 0.02 ]
        Minimum linewidth to plot to avoid large pdfs.

    color: string [ 'black' ]
        Color of lines.

    f(cweight) -> linewidth:
        (cmin, cmax) -> (lmax, lmin) if cmin < cweight < cmax

    """
    from matplotlib.collections import LineCollection
    mesh.createNeighbourInfos()

    lines = []
    linewidths = []

    bi = 0
    if cell_indices is None:
        # case 1/2: no node indices means we try plotting in implicit order of
        # boundaries in mesh
        boundary_ids = np.arange(len(cweight))

        bi = 0
        for bound in mesh.boundaries():
            if bound.rightCell() is not None:
                pos0 = np.array(bound.node(0).pos())[:2].tolist()
                pos1 = np.array(bound.node(1).pos())[:2].tolist()
                linew = lineWidthFromC(cweight[bi], cmin=cmin, cmax=cmax,
                                       lmin=lmin, lmax=lmax)

                if linew >= min_plot:
                    lines.append([pos0, pos1])
                    linewidths.append(linew)

                bi += 1

    else:
        # case 2/2: node indices are given, which means we plot each entry of
        # cweight between the nodes with the ids given in node_indices.
        assert len(cweight) == len(cell_indices[0])
        assert len(cweight) == len(cell_indices[1])
        cells = mesh.cells()
        boundary_ids = []
        for ci in range(len(cweight)):
            cell0 = cells[cell_indices[0][ci]]
            cell1 = cells[cell_indices[1][ci]]
            try:
                boundary_ids.append(
                    np.intersect1d([cell0.boundary(0).id(),
                                    cell0.boundary(1).id(),
                                    cell0.boundary(2).id()],
                                   [cell1.boundary(0).id(),
                                    cell1.boundary(1).id(),
                                    cell1.boundary(2).id()]
                                   )[0])
            except IndexError:
                print([cell0.boundary(0).id(),
                       cell0.boundary(1).id(),
                       cell0.boundary(2).id()],
                      [cell1.boundary(0).id(),
                       cell1.boundary(1).id(),
                       cell1.boundary(2).id()])

        print(boundary_ids)

        boundaries = mesh.boundaries()

        for ci, bi in enumerate(boundary_ids):
            bound = boundaries[bi]
            if bound.rightCell() is not None:
                pos0 = np.array(bound.node(0).pos())[:2].tolist()
                pos1 = np.array(bound.node(1).pos())[:2].tolist()
                linew = lineWidthFromC(cweight[ci], cmin=cmin, cmax=cmax,
                                       lmin=lmin, lmax=lmax)
                if linew >= min_plot:
                    lines.append([pos0, pos1])
                    linewidths.append(linew)

    ax.add_collection(
        LineCollection(lines, linewidths=linewidths,
                       color=color))