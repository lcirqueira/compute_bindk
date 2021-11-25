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
PATHDIC = cdef.FLD_PATH


conc_arr = np.arange(1 , DIL_END + 1 , 1)

bindk = np.load("{0}/bindk/{1:.2f}mM/bindk.{1:.2f}mM.npy".format(COMP_PATH , DIL_REF))
#bindk = np.loadtxt("{0}/bindk/{1:.2f}mM/bindk.{1:.2f}mM.dat".format(COMP_PATH , conc))


mean_n = []


for conc in conc_arr:
    conclignum = np.ceil( (LIGNUM[REF_CONC] / REF_CONC) * conc)

    densityUV = conclignum / (SYS_VOL / 100000)

    occnum = bindk.shape[0]

    probs = np.zeros((occnum, 2))
    probs[:,0] = np.arange(occnum)

    for num in range(occnum):
        probs[num,1] = (densityUV ** num) * bindk[num]

    probs[:,1] /= probs[:,1].sum()

    np.savetxt("projprobs.{:.2f}mM.dat".format(conc) , probs)

    meanlig = np.sum(np.prod(probs, axis=1))
    mean_n.append(meanlig)


meanout = np.zeros((len(conc_arr) , 2))
meanout[:,0] = conc_arr
meanout[:,1] = mean_n
np.savetxt("mean_n_allconc.dat" , meanout , header="Concentration (mM)  <n>" , fmt="%.1f %.4f")
