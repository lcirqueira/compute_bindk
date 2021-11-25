import numpy as np
from scipy import special
import sys , os
sys.path.append("/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/compute_bindk")
import computebindk_definitions as cdef
sys.path.append("../../")
import folderdefs as folddefs

RT = 0.001985 * 300 #kcal/mol 
BETA = 1/RT

#STDCONC = 1/1660 ##1/A^3
STDCONC = 1000/1660 ##1/nm^3


CONC = folddefs.FOLDER_CONC

BINDWORK = cdef.CONC_BINDWORK

DENS_SCALE = cdef.DENS_SCALE
REF_CONC = cdef.REF_CONC

LIGNUM = cdef.LIGNUM
NMAX = cdef.NMAX

LIPNUM = cdef.LIPNUM

SITES = cdef.SITES    
PATHDIC = cdef.FLD_PATH

conc_arr = cdef.CONCLIST


def get_Kn(dens , vol , netW , n):
    narr = np.arange(n + 1)

    stepKn = np.ones(narr.shape[0])
    #stepKn[1:] = (volume / narr[1:]) * np.exp(-BETA * (dwp_alpha + res_unsolv - bindwork) )
    stepKn[1:] = (vol / narr[1:]) * np.exp(-BETA * netW) 

    Kn = stepKn.astype("longdouble").cumprod()
    return(stepKn , Kn)


def get_dG(dens , stepKn , n):
    narr = np.arange(n + 1)

    stepdG = np.zeros_like(stepKn)
    stepdG[1:] = (-1/BETA) * np.log(stepKn[1:] * STDCONC)

    dG = stepdG.cumsum() 
    return(stepdG , dG)


def get_unsolv_work(h_0 , W_0 , n , lignmax):
    narr = np.arange(n + 1)
    frac = (lignmax - narr) / ((lignmax - narr) + LIPNUM)
    neg_frac = frac < 0

    unsolv_work = W_0 - h_0 + (h_0 * np.square(1 - frac))
    unsolv_work[neg_frac] = W_0

    #return(unsolv_work)
    return(unsolv_work.mean())



volume = np.loadtxt("{}/MSMS_wtreptow_volume.dat".format(cdef.MSMS_SURF_PATH))
UV_volume = volume / 100000

#bindwdata = np.loadtxt("{}/avg_stepbindw.dat".format(cdef.BINDW_PATH))
#meanbindw = bindwdata[0]


unsolv_dG_par = np.loadtxt("{}/analysis/partition_coefficient/inv-datafit-partition.dat".format(cdef.MEMB_PATH))
h0 , W0 = -np.round(unsolv_dG_par[0] , 3) , np.round(unsolv_dG_par[3] , 2)


dwp_alpha = np.loadtxt("{}/allconc_alpha.dat".format(cdef.DWP_PATH))
alphadic = {data[0] : data[1] for data in dwp_alpha}
alphadic = {1 : 0.0,
            25: 0.0,
            50: 0.0,
            75: -0.001,
            100: -0.004,
            150: -0.008}



with open("compute_params.dat" , "w") as fout:
    fout.write("""General compute params:
Volume = {} A^3
h0 = {} kcal/mol
W*(x=0) = {} kcal/mol

""".format(volume , h0 , W0))
    for i , conc in enumerate(conc_arr):
        if not os.path.exists("{:.2f}mM".format(conc)):
            os.mkdir("{:.2f}mM".format(conc))

        conc_lignum = LIGNUM[conc]

        bindw = BINDWORK[i]

        conc_dwp_alpha = alphadic[conc]

        unsolv_work = get_unsolv_work(h0 , W0 , NMAX , conc_lignum)

        net_work = conc_dwp_alpha + unsolv_work - bindw


        stepbindk , bindk = get_Kn(conc , UV_volume , net_work , NMAX)
        np.savetxt("{0:.2f}mM/stepbindk.{0:.2f}mM.dat".format(conc) , stepbindk)
        np.savetxt("{0:.2f}mM/bindk.{0:.2f}mM.dat".format(conc) , bindk)
        np.save("{0:.2f}mM/bindk.{0:.2f}mM".format(conc) , bindk)


        stepdeltaG , deltaG = get_dG(conc , stepbindk , NMAX)
        #np.savetxt("{0:.2f}mM/stepdeltaG.{0:.2f}mM.dat".format(conc) , stepdeltaG)
        #np.savetxt("{0:.2f}mM/deltaG.{0:.2f}mM.dat".format(conc) , deltaG)


        fout.write("""{}mM compute bindk parameters:
W** = {:.4f} kcal/mol
alpha = {:.4f}
mean W* = {:.4f}
net W = {:.4f}

""".format(conc , bindw , conc_dwp_alpha , unsolv_work , net_work))
