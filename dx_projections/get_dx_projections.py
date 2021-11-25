import numpy as np
import gridData as gd
from scipy import special
import sys , os
sys.path.append("/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/compute_bindk")
import computebindk_definitions as cdef
sys.path.append("../../")
import folderdefs as folddefs
np.seterr(all='raise')

CONC = folddefs.FOLDER_CONC

REF_CONC = cdef.REF_CONC

DIL_REF = cdef.DIL_REF
DIL_END = cdef.DIL_END

COMP_PATH = cdef.COMPUTE_PATH

FLDSITE = "tmptn"    
PATHDIC = cdef.FLD_PATH

LIGNUM = cdef.LIGNUM

FLD_CONC = cdef.CONCLIST[1:]

#PROJ_CONC_ARR = [1,]
PROJ_CONC_ARR = np.arange(1 , DIL_END + 1 , 1)


dx_grids = {}
clusterprobs = {}

xlist = []
ylist = []
zlist = []

for occn in range(1 , LIGNUM[REF_CONC] + 1):
    dx_grids[occn] = []
    clusterprobs[occn] = []

    for fldconc in FLD_CONC:
        path = PATHDIC[fldconc]

        clprobs = np.loadtxt("{0}/analysis/allptn/clusters/probs.cluster.tmptn.dat".format(path))

        dxpath = "{0}/analysis/allptn/clusters/volmap.cluster.{1}.{2}.dx".format(path , FLDSITE , occn)
        if os.path.exists(dxpath):
            clusterprobs[occn].append(clprobs[occn - 1, 1])

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
    occ_grid_probs = clusterprobs[occn]
    occ_probs_sum = np.sum(occ_grid_probs)

    occ_grid_list = dx_grids[occn]

    occgrid = dxmap.resample(edgeslist)
    occgrid.grid = np.zeros_like(occgrid.grid)

    for tmpgrid , tmpprobs in zip(occ_grid_list , occ_grid_probs):
        if occ_probs_sum != 0:
            occgrid.grid += tmpgrid.resample(occgrid).grid * (tmpprobs / occ_probs_sum)

    allndic[occn] = occgrid

    if occgrid.grid.sum() != 0:
        occgrid.export("occmean.{}.dx".format(occn))


conc_mean_num = np.loadtxt("{}/dilute_projections/mean_n_allconc.dat".format(COMP_PATH))
mean_num_dic = {dens : avgn for dens , avgn in conc_mean_num}


with open("mean_n_projections.dat" , "w") as fout:
    fout.write("#Conc   <n>\n")
    for c , projconc in enumerate(PROJ_CONC_ARR):
        probs = np.loadtxt("{}/dilute_cluster/probs.cluster.{}.{:.2f}mM.dat".format(COMP_PATH , FLDSITE , projconc))

        occ_states = probs[:,0]

        projgrid = occgrid.resample(occgrid)
        projgrid.grid = np.zeros_like(projgrid.grid)

        for occ in occ_states:
            occmap = allndic[occ]
            mean_n = mean_num_dic[projconc]

            projgrid.grid += (mean_n * occmap.resample(projgrid).grid * probs[int(occ) - 1 , 1])


        projgrid.export("proj.{}.{}mM.dx".format(FLDSITE , projconc))
        fout.write("{:.0f}   {:.4f}\n".format(projconc , projgrid.grid.sum() ))
