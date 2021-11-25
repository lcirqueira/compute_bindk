import numpy as np

REF_CONC = 150 #mM
#CONC_ARR = range(1 , 151 , 0.1)

BINDWORK = 4.57
CONC_BINDWORK = [5.02 , 5.02 , 5.00 , 4.9 , 4.8 , 4.64]

NET_WORK = -0.77

LIGNUM = {150 : 174,
          100 : 116,
           75 : 87,
           50 : 58,
           25 : 29,
           1  : 2} 
#LIGNUM = {150 : 72,
          #100 : 72,
           #75 : 72,
           #50 : 72,
           #25 : 72,
           #1  : 72}


FLD_PATH = {150 : "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane",
            100 : "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevo100mM",
             75 : "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevo75mM",
             50 : "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevo50mM",
             25  : "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevo25mM",
             1  : "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevo1mM",
             0  : "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/jhosoume",}

SITES = ("tmptn" , "centralpore" , "wt")    

COMPUTE_PATH = "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/compute_bindk"

BINDW_PATH = "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/lstock/bindw"

PTNSURF_PATH = "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/protein_surface/surf_volume"
MSMS_SURF_PATH = "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/wtreptowMSMS/surfvolume"

MEMB_PATH = "/home/lcirqueira/Simulations/membrane/flooding/sevoflurane"
FRACTION =  [ 0.0043 , 0.06 , 0.11 , 0.16 , 0.21 , 0.29]
CONCLIST = [1 , 25 , 50 , 75 , 100 , 150]
CONCCOLOR = {1 : "#009f2b" , 25 : "#01a9f5" , 50 : "#00abac" , 75 : "#a8a9ab" , 100 : "#feb3a9" , 150 : "#e4201f"}

LIPNUM = 426
SYSVOL = 1236104

DWP_PATH = "/home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/ptn_struct/deltaWp"

DENS_SCALE = 1.083 #from /home/lcirqueira/Simulations/ionchannel/kv1/kv1.12/flooding/sevoflurane/analysis/wtreptowMSMS/res_volume
