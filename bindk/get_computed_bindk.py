import numpy as np
from scipy import special
import sys , os
sys.path.append("/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/compute_bindk")
import computebindk_definitions as cdef
sys.path.append("../../")
import folderdefs as folddefs

RT = 0.001985 * 300 #kcal/mol 
BETA = 1/RT

STDCONC = 1/1660

CONC = folddefs.FOLDER_CONC

BINDWORK = cdef.BINDWORK

DENS_SCALE = cdef.DENS_SCALE
REF_CONC = cdef.REF_CONC

LIGNUM = cdef.LIGNUM

SITES = cdef.SITES    
PATHDIC = cdef.FLD_PATH

conc_arr = cdef.CONCLIST


def get_Kn(dens , n):
    narr = np.arange(n + 1)
    dwp_alpha = alphadic[dens]
    res_unsolv = unsolvdic[dens] 

    stepKn = np.ones(narr.shape[0])
    stepKn[1:] = (volume / narr[1:]) * np.exp(-BETA * (dwp_alpha + res_unsolv - meanbindw) )

    Kn = stepKn.astype("longdouble").cumprod()
    return(stepKn , Kn)


def get_dG(dens , n):
    stepKn , Kn = get_Kn(dens , n)
    narr = np.arange(n + 1)

    stepdG = np.zeros_like(stepKn)
    stepdG[1:] = (-1/BETA) * np.log(stepKn[1:] * STDCONC)

    dG = stepdG.cumsum() 
    return(stepdG , dG)



volume = np.loadtxt("{}/MSMS_wtreptow_volume.dat".format(cdef.MSMS_SURF_PATH))

#bindwdata = np.loadtxt("{}/avg_stepbindw.dat".format(cdef.BINDW_PATH))
#meanbindw = bindwdata[0]
meanbindw = BINDWORK


unsolv_dG = np.loadtxt("{}/analysis/partition_coefficient/fep_energies.dat".format(cdef.MEMB_PATH))
unsolvdic = {data[0] : data[1] for data in unsolv_dG}


dwp_alpha = np.loadtxt("{}/allconc_alpha.dat".format(cdef.DWP_PATH))
alphadic = {data[0] : data[1] for data in dwp_alpha}

for conc in conc_arr:
    if not os.path.exists("{:.2f}mM".format(conc)):
        os.mkdir("{:.2f}mM".format(conc))

    stepbindk , bindk = get_Kn(conc , LIGNUM[conc])
    np.savetxt("{0:.2f}mM/stepbindk.{0:.2f}mM.dat".format(conc) , stepbindk)
    np.savetxt("{0:.2f}mM/bindk.{0:.2f}mM.dat".format(conc) , bindk)
    np.save("{0:.2f}mM/bindk.{0:.2f}mM".format(conc) , bindk)

    stepdeltaG , deltaG = get_dG(conc , LIGNUM[conc])
    np.savetxt("{0:.2f}mM/stepdeltaG.{0:.2f}mM.dat".format(conc) , stepdeltaG)
    np.savetxt("{0:.2f}mM/deltaG.{0:.2f}mM.dat".format(conc) , deltaG)


with open("compute_params.dat" , "w") as fout:
    fout.write("""Compute bindk parameters:
W** = {:.3f} kcal/mol

Volume = {} A^3

W* = {}

alpha = {}
""".format(BINDWORK , volume , unsolvdic , alphadic))
