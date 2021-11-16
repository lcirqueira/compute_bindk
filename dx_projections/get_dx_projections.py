import numpy as np
import gridData as gd
from scipy import special
import sys , os
sys.path.append("/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/compute_bindk")
import computebindk_definitions as cdef
sys.path.append("../../")
import folderdefs as folddefs

CONC = folddefs.FOLDER_CONC

REF_CONC = cdef.REF_CONC

COMP_PATH = cdef.COMPUTE_PATH

FLDSITE = "tmptn"    
PATHDIC = cdef.FLD_PATH

LIGNUM = cdef.LIGNUM

FLD_CONC = cdef.CONCLIST[1:]

PROJ_CONC_ARR = cdef.CONCLIST
PROJ_CONC_ARR = [1,]


dx_grids = {}

xlist = []
ylist = []
zlist = []

for occn in range(1 , LIGNUM[REF_CONC] + 1):
    dx_grids[occn] = []

    for fldconc in FLD_CONC:
        path = PATHDIC[fldconc]

        dxpath = "{0}/analysis/allptn/clusters/volmap.cluster.{1}.{2}.dx".format(path , FLDSITE , occn)
        if os.path.exists(dxpath):
            dxmap = gd.Grid(dxpath)
            dx_grids[occn].append(dxmap)

            xlist += [dxmap.edges[0][0] , dxmap.edges[0][-1]]
            ylist += [dxmap.edges[1][0] , dxmap.edges[1][-1]]
            zlist += [dxmap.edges[2][0] , dxmap.edges[2][-1]]



xedges = min(xlist) , max(xlist)
yedges = min(ylist) , max(ylist)
zedges = min(zlist) , max(zlist)

edgeslist = [ np.arange(xedges[0] , xedges[1] + 1) , np.arange(yedges[0] , yedges[1] + 1) , np.arange(zedges[0] , zedges[1] + 1) ]

allndic = {}

for occn in range(1 , LIGNUM[REF_CONC] + 1):
    occ_grid_list = dx_grids[occn]

    occgrid = dxmap.resample(edgeslist)
    occgrid.grid = np.zeros_like(occgrid.grid)

    for tmpgrid in occ_grid_list:
        occgrid.grid += tmpgrid.resample(occgrid).grid * (1/len(occ_grid_list))

    allndic[occn] = occgrid

    occgrid.export("occmean.{}.dx".format(occn))

with open("mean_n_projections.dat" , "w") as fout:
    fout.write("#Conc   <n>\n")
    for projconc in PROJ_CONC_ARR:
        probs = np.loadtxt("{0}/projections/projprobs.{1:.2f}mM.dat".format(COMP_PATH , projconc))

        occ_states = probs[1:,0]

        projgrid = occgrid.resample(occgrid)
        projgrid.grid = np.zeros_like(projgrid.grid)

        for occ in occ_states:
            occmap = allndic[occ]

            projgrid.grid += occmap.resample(projgrid).grid * probs[int(occ) , 1]


        projgrid.export("proj.{}.{}mM.dx".format(FLDSITE , projconc))
        fout.write("{:.0f}   {:.4f}\n".format(projconc , projgrid.grid.sum() ))
