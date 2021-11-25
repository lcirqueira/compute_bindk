import numpy as np
from scipy import special
import sys , os
sys.path.append("/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/compute_bindk")
import computebindk_definitions as cdef
sys.path.append("../../")
import folderdefs as folddefs

RT = 0.001985 * 300 #kcal/mol 
BETA = 1/RT

CONC = folddefs.FOLDER_CONC

DENS_SCALE = cdef.DENS_SCALE
REF_CONC = cdef.REF_CONC

DIL_REF = cdef.DIL_REF
DIL_END = cdef.DIL_END

LIGNUM = cdef.LIGNUM
SYS_VOL = cdef.SYS_VOL

COMP_PATH = cdef.COMPUTE_PATH

SITES = cdef.SITES    
SITE = "tmptn"
PATHDIC = cdef.FLD_PATH


conc_arr = np.arange(0.5 , DIL_END + 0.5 , 0.5)


dil_ref_probs = np.loadtxt("{}/analysis/allptn/clusters/probs.cluster.tmptn.dat".format(PATHDIC[DIL_REF]))
dil_end_probs = np.loadtxt("{}/analysis/allptn/clusters/probs.cluster.tmptn.dat".format(PATHDIC[DIL_END]))


delta_probs = (dil_end_probs[:,1] - dil_ref_probs[:,1]) / (DIL_END - DIL_REF)


for conc in conc_arr:
    if conc < DIL_REF:
        conc_cluster_probs = dil_ref_probs.copy()
    else:
        conc_cluster_probs = dil_ref_probs.copy()
        conc_cluster_probs[:,1] = (delta_probs * (conc - DIL_REF)) + dil_ref_probs[:,1]

    np.savetxt("probs.cluster.{}.{:.2f}mM.dat".format(SITE , conc) , conc_cluster_probs)
