####### This python scipt helps analyze the trajectories from Cytosim

######## Imports ########
import math
import numpy as np
import os
import pandas as pd
import icosphere as ic
from datetime import datetime
from collections import Counter
from collections import defaultdict
import os
import subprocess

class replicate:
    ridx = 0;
    Nsnaps = 0;
    snap = [];
    timevector = []
    def __init__(self, rid):
        self.ridx = rid
        self.Nsnaps = 0
        self.snap = [];
        self.timevector =[];
    def addsnapshot(self):
        self.snap.append(snapshot(1))

class snapshot:
    sidx = 0;
    filcoord=[]
    crossboundcoord=[];
    crossboundblobid = [];
    handid = [];
    dimerid = [];
    boundid = [];
    partnerid = [];
    crossboundfilid = [];
    solidcoord = [];
    selfblobid = []
    def __init__(self, sid):
        self.sidx = sid
        self.filcoord = [];
        self.crossboundcoord = [];
        self.crossboundblobid = [];
        self.handid = [];
        self.dimerid = [];
        self.boundid = [];
        self.partnerid = [];
        self.crossboundfilid = [];
        self.solidcoord = [];
        self.selfblobid = [];

def readcytosimtraj (*args):
    tag = ''
    filepath = args[0]
    filename = 'objects.cmo'
    if(len(args)==2):
        tag= args[1]
    
    # Get number of solids
    test = subprocess.run(["bash", "./get_nsolids.sh", "-f", filepath], capture_output=True)
    nsolids = int(test.stdout)
    print("Nsolids in trajectory = "+str(nsolids), flush=True)
    r=[];
    Nruns =1 
    for  i in range(0,Nruns):
        r.append(replicate(i))
    printstatus = False; 
    recordStatus = False; 
    ridx = 0;
    sidx=-1;
    # Open file
    fptr = open(filepath+filename,'r', encoding="utf8", errors="ignore")
    for line in fptr:
        if(printstatus):
            print(line)

        if('time' in line):
            t = float((line.split(' '))[1])
            r[ridx].timevector.append(t)
            r[ridx].addsnapshot();
            sidx = sidx + 1
            r[ridx].snap[sidx].solidcoord = [None]*nsolids
        
        if('#section fiber' in line): 
            fcoord=[];
            line = fptr.readline()
            if(printstatus):
                print(line)
            while(not('section' in line)):
                if(printstatus):
                    print(line)
                if('f2' in line or 'f 2' in line):
                    recordStatus = False;
                if('f1' in line or 'f 1' in line):
                    if(fcoord):
                        r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                        fcoord=[];
                    recordStatus = True;
                elif(' ' in line[0] and recordStatus):
                    line = line.strip()
                    cstring = line.split(' ')
                    fcoord.append([float(cstring[0]), float(cstring[1]), float(cstring[2])])
                line = fptr.readline()
            #when it exists, if fcoord has not been recorded, record it.
            if(fcoord):
                r[ridx].snap[sidx].filcoord.append(np.array(fcoord))
                fcoord=[];
        
        # crosslinker hands
        if('#section solid' in line and 'filonly' != tag):
            line = fptr.readline()
            line = line.strip()
            if(len(line) == 0):
                line = fptr.readline()
            while(not('section' in line)):
                if(printstatus):
                    print(line)
                #Get solid ID
                line = line.strip()
                cstring = line.split(' ')
                if(len(cstring) == 2):
                    cstring = cstring[0].split(':') 
                else:
                    cstring = cstring[1].split(':') 
                solidid = int(cstring[1]) 
                solidcoord = [];
                line = fptr.readline()
                while(line[0]!='d' and not('section' in line)):
                    #Get coord
                    line = line.strip()
                    cstring = line.split(' ')
                    solidcoord.append([float(cstring[0]), float(cstring[1]), float(cstring[2])])
                    line = fptr.readline()
                if(solidcoord):
                    r[ridx].snap[sidx].solidcoord[solidid-1] = np.array(solidcoord[0][:])
            
            # Bound crosslinker hands
            cstring_temp = []
            if('#section single A' in line and 'filonly' != tag):
                line = fptr.readline()
                if(printstatus):
                    print(line)
                while(not('section' in line)):
                    if(printstatus):
                        print(line)
                    if line[0]=='w':
                        line = line.strip()
                        cstring = line.split(' ')
                        if(len(cstring) <= 8):
                            cstring_temp2 = cstring[0].split(':')
                            cstring_temp = float(cstring_temp2[1])
                            if(len(cstring) == 7):
                                r[ridx].snap[sidx].partnerid.append(int(cstring[1][1:len(cstring[1])]))
                                r[ridx].snap[sidx].selfblobid.append(int(cstring[2][1:len(cstring[2])]))
                                r[ridx].snap[sidx].dimerid.append(int(cstring_temp))
                            elif(len(cstring) == 8):
                                r[ridx].snap[sidx].crossboundfilid.append(int(cstring[1][1:len(cstring[1])]))
                                r[ridx].snap[sidx].crossboundblobid.append(int(cstring[3][1:len(cstring[3])]))
                                r[ridx].snap[sidx].boundid.append(int(cstring_temp))
                            else:
                                r[ridx].snap[sidx].crossboundfilid.append(int(cstring[2][1:len(cstring[2])]))
                                r[ridx].snap[sidx].crossboundblobid.append(int(cstring[4][1:len(cstring[4])]))
                        elif(len(cstring) > 8):
                            cstring_temp2 = cstring[1].split(':')
                            cstring_temp = float(cstring_temp2[1])
                            if(len(cstring) == 10):
                                r[ridx].snap[sidx].partnerid.append(int(cstring[3]))
                                r[ridx].snap[sidx].selfblobid.append(int(cstring[5]))
                                r[ridx].snap[sidx].dimerid.append(int(cstring_temp))
                            elif(len(cstring) == 11):
                                r[ridx].snap[sidx].crossboundfilid.append(int(cstring[3]))
                                r[ridx].snap[sidx].crossboundblobid.append(int(cstring[6]))
                                r[ridx].snap[sidx].boundid.append(int(cstring_temp))
                    line = fptr.readline()
    fptr.close()
    print('Number of snapshots='+str(len(r[ridx].snap)))
    fptr.close()
    return r

def getfracoccupiedtimeseries(N, dirname, meshvar, Nlists, unitnormal, Ntriangles, r, deltasnap, outfilename):
    Nsnaps = len(r[0].snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    fvec = np.zeros((2,Ndatapoints))
    fvec[0][:]=np.arange(0,Nsnaps,deltasnap)
    for SREF in range(0,Nsnaps,deltasnap):
        print(SREF,flush=True)
        filcoord = r[0].snap[SREF].filcoord 
        actincounter = ic.generatedensityfield(meshvar, unitnormal, filcoord, False, ic.SearchAlgoType.LISTEDSEARCH, Nlists)
        #Order Parameter #3
        is_all_zero = np.argwhere((actincounter == 0.0))
        fvec[1][SREF] = 1-len(is_all_zero)/Ntriangles
    pd.DataFrame(fvec, index=['Time','Frac_occupied']).to_csv(outfilename)

def getvaspprops(rset, deltasnap, N, Rval, outputfile):
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    bar_min_snap = int(0.95*Nsnaps)
    bar_data_raw = np.array([])
    collectstatus = False
    datawritestatus = False
    datamatrix = np.zeros((19,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    crosslinker_n_fil_bar = [[],[],[],[]]

    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        nvasp_tmp = np.array([])
        nvalency_tmp = np.array([])
        Nfil_per_solid_tmp = np.array([])
        if SREF>=bar_min_snap:
            collectstatus = True
        countermat = np.zeros((4,len(rset)))
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                # The two variables below together represent information of each bond in the network.
                crossboundfilid = np.array(rset[runidx].snap[SREF].crossboundfilid)
                crossboundblobid = rset[runidx].snap[SREF].crossboundblobid
                if(len(crossboundblobid)):
                    datawritestatus = True
                unique, counts = np.unique(crossboundblobid, return_counts=True)
                nvasp_tmp = np.append(nvasp_tmp,len(unique))
                nvalency_tmp = np.append(nvalency_tmp,np.array(counts))
                # Get the filament IDs for each solid and calculate Nhands/Nfil
                # This helps you understand the number of hands in each solid that are bound to the same filament.
                # The distribution of Nhands/Nfil - closer to 1 suggests that each bound hand is bound to a different filament
                # >1 suggests that multiple hands of the solid are bound to the same fil ID
                # <1 is not possible.
                for ublobid in enumerate(unique):
                    #Find the locs using argwhere
                    locs = np.argwhere(crossboundblobid==ublobid)
                    #Corresponding filament IDs that share this 
                    #Crosslinker in the bond
                    filid_solid = crossboundfilid[locs]
                    #Get filament IDs and the counts
                    unique_fid, counts_fil = np.unique(filid_solid, return_counts=True)
                    #Get the total number of unique filaments
                    lval = len(unique_fid)
                    countermat[lval-1][runidx] = countermat[lval-1][runidx] + 1
                    Nfil_per_solid_tmp = np.append(Nfil_per_solid_tmp,lval)
                if collectstatus:
                    bar_data_raw = np.append(bar_data_raw, Nfil_per_solid_tmp)
                    for pos in range(0,4):
                        crosslinker_n_fil_bar[pos].append(countermat[pos][runidx])

        meanvec = np.mean(countermat,1)
        sstdvec = np.std(countermat,1)

        if(datawritestatus):
            datamatrix[1][scounter]= np.mean(nvasp_tmp)
            datamatrix[2][scounter]= np.std(nvasp_tmp)
            datamatrix[3][scounter]= np.mean(nvalency_tmp)
            datamatrix[4][scounter]= np.std(nvalency_tmp)
            datamatrix[5][scounter]= np.mean(Nfil_per_solid_tmp)
            datamatrix[6][scounter]= np.std(Nfil_per_solid_tmp)
            datamatrix[9][scounter]= meanvec[0]
            datamatrix[10][scounter]= sstdvec[0]
            datamatrix[11][scounter]= meanvec[1]
            datamatrix[12][scounter]= sstdvec[1]
            datamatrix[13][scounter]= meanvec[2]
            datamatrix[14][scounter]= sstdvec[2]
            datamatrix[15][scounter]= meanvec[3]
            datamatrix[16][scounter]= sstdvec[3]

        scounter = scounter + 1 
        if(len(bar_data_raw)):
            datamatrix[7][0] = np.mean(bar_data_raw)
            datamatrix[8][0] =  np.std(bar_data_raw)
            for pos in range(0,4):
                datamatrix[17][pos] = np.mean(crosslinker_n_fil_bar[pos])
                datamatrix[18][pos] =  np.std(crosslinker_n_fil_bar[pos])
        
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    #Nfpc = Number of filaments per crosslinker
    pd.DataFrame(datamatrix, index=['Time', 'Ncross_mean','Ncross_std','Valency_mean','Valency_std',
                                    'Nfil_per_cross_mean','Nfil_per_cross_std',
                                    'Bar_Nfil_per_cross_mean','Bar_Nfil_per_cross_std',
                                    'Mean_Nfpc_1_fil','Std_Npc_1_fil',
                                    'Mean_Nfpc_2_fil','Std_Nfpc_2_fil',
                                    'Mean_Nfpc_3_fil','Std_Nfpc_3_fil',
                                    'Mean_Nfpc_4_fil','Std_Nfpc_4_fil',
                                    'Bar_Mean_Nfpc','Bar_Std_Nfpc'
                                    ]).to_csv(outputfile)

def getvasppropsLpdVASP(rset, deltasnap, N, Rval, outputfile, ratio1):
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    bar_min_snap = int(0.95*Nsnaps)
    bar_data_bivalent_raw = np.array([])
    bar_data_tetravalent_raw = np.array([])
    collectstatus = False
    datawritestatus = False
    datamatrix = np.zeros((33,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    crosslinker_n_fil_bar_bivalent = [[],[],[],[]]
    crosslinker_n_fil_bar_tetravalent = [[],[],[],[]]

    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        nvasp_tmp = np.array([])
        nvalency_tmp = np.array([])
        Nfil_per_solid_bivalent_tmp = np.array([])
        Nfil_per_solid_tetravalent_tmp = np.array([])
        if SREF>=bar_min_snap:
            collectstatus = True
        countermat_bivalent = np.zeros((4,len(rset)))
        countermat_tetravalent = np.zeros((4,len(rset)))
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                # The two variables below together represent information of each bond in the network.
                crossboundfilid = np.array(rset[runidx].snap[SREF].crossboundfilid)
                crossboundblobid = rset[runidx].snap[SREF].crossboundblobid
                if(len(crossboundblobid)):
                    datawritestatus = True
                unique, counts = np.unique(crossboundblobid, return_counts=True)
                nvasp_tmp = np.append(nvasp_tmp,len(unique))
                nvalency_tmp = np.append(nvalency_tmp,np.array(counts))
                # Get the filament IDs for each solid and calculate Nhands/Nfil
                # This helps you understand the number of hands in each solid that are bound to the same filament.
                # The distribution of Nhands/Nfil - closer to 1 suggests that each bound hand is bound to a different filament
                # >1 suggests that multiple hands of the solid are bound to the same fil ID
                # <1 is not possible.
                for ublobid in enumerate(unique):
                    if ublobid <= int(ratio1):
                        #Find the locs using argwhere
                        locs = np.argwhere(crossboundblobid==ublobid)
                        #Corresponding filament IDs that share this crosslinker in the bond
                        filid_solid = crossboundfilid[locs]
                        #Get filament IDs and the counts
                        unique_fid = np.unique(filid_solid, return_counts=True)
                        #Get the total number of unique filaments
                        lval = len(unique_fid)
                        countermat_bivalent[lval-1][runidx] = countermat_bivalent[lval-1][runidx] + 1
                        Nfil_per_solid_bivalent_tmp = np.append(Nfil_per_solid_bivalent_tmp,lval)
                        
                    else:
                        #Find the locs using argwhere
                        locs = np.argwhere(crossboundblobid==ublobid)
                        #Corresponding filament IDs that share this 
                        #Crosslinker in the bond
                        filid_solid = crossboundfilid[locs]
                        #Get filament IDs and the counts
                        unique_fid = np.unique(filid_solid, return_counts=True)
                        #Get the total number of unique filaments
                        lval = len(unique_fid)
                        countermat_tetravalent[lval-1][runidx] = countermat_tetravalent[lval-1][runidx] + 1
                        Nfil_per_solid_tetravalent_tmp = np.append(Nfil_per_solid_tetravalent_tmp,lval)

                if collectstatus:
                    bar_data_bivalent_raw = np.append(bar_data_bivalent_raw, Nfil_per_solid_bivalent_tmp)
                    bar_data_tetravalent_raw = np.append(bar_data_tetravalent_raw, Nfil_per_solid_tetravalent_tmp)
                    for pos in range(0,4):
                        crosslinker_n_fil_bar_bivalent[pos].append(countermat_bivalent[pos][runidx])
                        crosslinker_n_fil_bar_tetravalent[pos].append(countermat_tetravalent[pos][runidx])

        meanvec_bivalent = np.mean(countermat_bivalent,1)
        sstdvec_bivalent = np.std(countermat_bivalent,1)
        meanvec_tetravalent = np.mean(countermat_tetravalent,1)
        sstdvec_tetravalent = np.std(countermat_tetravalent,1)

        if(datawritestatus):
            datamatrix[1][scounter]= np.mean(nvasp_tmp)
            datamatrix[2][scounter]= np.std(nvasp_tmp)
            datamatrix[3][scounter]= np.mean(nvalency_tmp)
            datamatrix[4][scounter]= np.std(nvalency_tmp)
            datamatrix[5][scounter]= np.mean(Nfil_per_solid_bivalent_tmp)
            datamatrix[6][scounter]= np.std(Nfil_per_solid_bivalent_tmp)
            datamatrix[9][scounter]= meanvec_bivalent[0]
            datamatrix[10][scounter]= sstdvec_bivalent[0]
            datamatrix[11][scounter]= meanvec_bivalent[1]
            datamatrix[12][scounter]= sstdvec_bivalent[1]
            datamatrix[13][scounter]= meanvec_bivalent[2]
            datamatrix[14][scounter]= sstdvec_bivalent[2]
            datamatrix[15][scounter]= meanvec_bivalent[3]
            datamatrix[16][scounter]= sstdvec_bivalent[3]
            
            datamatrix[19][scounter]= np.mean(Nfil_per_solid_tetravalent_tmp)
            datamatrix[20][scounter]= np.std(Nfil_per_solid_tetravalent_tmp)
            datamatrix[23][scounter]= meanvec_tetravalent[0]
            datamatrix[24][scounter]= sstdvec_tetravalent[0]
            datamatrix[25][scounter]= meanvec_tetravalent[1]
            datamatrix[26][scounter]= sstdvec_tetravalent[1]
            datamatrix[27][scounter]= meanvec_tetravalent[2]
            datamatrix[28][scounter]= sstdvec_tetravalent[2]
            datamatrix[29][scounter]= meanvec_tetravalent[3]
            datamatrix[30][scounter]= sstdvec_tetravalent[3]

        scounter = scounter + 1 
        if(len(bar_data_bivalent_raw)):
            datamatrix[7][0] = np.mean(bar_data_bivalent_raw)
            datamatrix[8][0] =  np.std(bar_data_bivalent_raw)
            for pos in range(0,4):
                datamatrix[17][pos] = np.mean(crosslinker_n_fil_bar_bivalent[pos])
                datamatrix[18][pos] =  np.std(crosslinker_n_fil_bar_bivalent[pos])
        if(len(bar_data_tetravalent_raw)):
            datamatrix[21][0] = np.mean(bar_data_tetravalent_raw)
            datamatrix[22][0] =  np.std(bar_data_tetravalent_raw)
            for pos in range(0,4):
                datamatrix[31][pos] = np.mean(crosslinker_n_fil_bar_tetravalent[pos])
                datamatrix[32][pos] =  np.std(crosslinker_n_fil_bar_tetravalent[pos])
        
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    #Nfpc = Number of filaments per crosslinker
    pd.DataFrame(datamatrix, index=['Time', 'Ncross_mean','Ncross_std','Valency_mean','Valency_std',
                                    'Nfil_per_bivalent_cross_mean','Nfil_per_bivalent_cross_std',
                                    'Bar_Nfil_per_bivalent_cross_mean','Bar_Nfil_per_bivalent_cross_std',
                                    'Mean_Nfpbc_1_fil','Std_Nfpbc_1_fil',
                                    'Mean_Nfpbc_2_fil','Std_Nfpbc_2_fil',
                                    'Mean_Nfpbc_3_fil','Std_Nfpbc_3_fil',
                                    'Mean_Nfpbc_4_fil','Std_Nfpbc_4_fil',
                                    'Bar_Mean_Nfpbc','Bar_Std_Nfpbc',
                                    'Nfil_per_tetravalent_cross_mean','Nfil_per_tetravalent_cross_std',
                                    'Bar_Nfil_per_tetravalent_cross_mean','Bar_Nfil_per_tetravalent_cross_std',
                                    'Mean_Nfptc_1_fil','Std_Nfptc_1_fil',
                                    'Mean_Nfptc_2_fil','Std_Nfptc_2_fil',
                                    'Mean_Nfptc_3_fil','Std_Nfptc_3_fil',
                                    'Mean_Nfptc_4_fil','Std_Nfptc_4_fil',
                                    'Bar_Mean_Nfptc','Bar_Std_Nfptc'
                                    ]).to_csv(outputfile)

def getDynamicDimerprops(rset, deltasnap, N, Rval, outputfile):
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    datawritestatus = False
    datamatrix = np.zeros((10,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)

    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)
        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                # The two variables below together represent information of each bond in the network.
                # crossboundfilid = np.array(rset[runidx].snap[SREF].crossboundfilid)
                crossboundblobid = rset[runidx].snap[SREF].crossboundblobid
                partnerid = rset[runidx].snap[SREF].partnerid
                # selfblobid = rset[runidx].snap[SREF].selfblobid
                
                if(len(crossboundblobid)):
                    datawritestatus = True
                dimerid = rset[runidx].snap[SREF].dimerid
                partneridblob = rset[runidx].snap[SREF].partnerid
                partnerid = [2*(number)-1 for number in partneridblob]
                boundid = rset[runidx].snap[SREF].boundid

                partnermap = np.column_stack((dimerid, partnerid))
                partnermap_dict = dict(partnermap)

                nHands = 2*len(rset[runidx].snap[SREF].solidcoord)

                all_dimer_hands = [number for number in range(1,nHands+1) if number % 2 != 0]
                monomerid = [number for number in all_dimer_hands if number not in dimerid]

                all_actin_hands = [number for number in range(2,nHands+2) if number % 2 == 0]
                freeid = [number for number in all_actin_hands if number not in boundid]

                corresponding_boundid = [number-1 for number in boundid]
                bound_dimerid = [number for number in dimerid if number in corresponding_boundid]
                bound_monomerid = [number for number in monomerid if number in corresponding_boundid]
                bound_partnerid = [number for number in partnerid if number in corresponding_boundid]

                corresponding_freeid = [number-1 for number in freeid]
                free_dimerid = [number for number in dimerid if number in corresponding_freeid]
                free_monomerid = [number for number in monomerid if number in corresponding_freeid]
                free_partnerid = [number for number in partnerid if number in corresponding_freeid]
                
                bound_dimer_remapped = [partnermap_dict.get(n, n) for n in bound_dimerid]
                bd_bp = [number for number in bound_dimer_remapped if number in bound_partnerid]
                bd_fp = [number for number in bound_dimer_remapped if number in free_partnerid]

                free_dimer_remapped = [partnermap_dict.get(n, n) for n in free_dimerid]
                fd_bp = [number for number in free_dimer_remapped if number in bound_partnerid]
                fd_fp = [number for number in free_dimer_remapped if number in free_partnerid]

                # Valency Lists
                Monomer_0 = free_monomerid
                Monomer_1 = bound_monomerid
                Dimer_0 = fd_fp
                Dimer_1 = bd_fp + fd_bp
                Dimer_2 = bd_bp

        if(datawritestatus):
            # Not looking across the dimer bond
            datamatrix[1][scounter]= len(bound_dimerid)
            datamatrix[2][scounter]= len(bound_monomerid)
            datamatrix[3][scounter]= len(free_dimerid)
            datamatrix[4][scounter]= len(free_monomerid)
            # Looking across the dimer bond
            datamatrix[5][scounter]= len(Monomer_0)
            datamatrix[6][scounter]= len(Monomer_1)
            datamatrix[7][scounter]= len(Dimer_0)
            datamatrix[8][scounter]= len(Dimer_1)
            datamatrix[9][scounter]= len(Dimer_2)

        scounter = scounter + 1 
        
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    pd.DataFrame(datamatrix, index=['Time', 'Bound_Dimer', 'Bound_Monomer', 'Free_Dimer', 'Free_Monomer', 
                                    'Monomer_0', 'Monomer_1', 'Dimer_0', 'Dimer_1', 'Dimer_2']).to_csv(outputfile)

def find_chains_rings_and_doubles(mapping_dict):
    visited = set()
    chains = []
    rings = []
    doubles = set()

    def trace_chain(start_id):
        stack = [(start_id, None)]
        path = []
        ring_detected = False
        ring_start = None
        while stack:
            current_id, prev_id = stack.pop()
            if current_id in path:
                ring_start = path.index(current_id)
                ring_detected = True
                break
            if current_id in visited:
                continue
            visited.add(current_id)
            path.append(current_id)

            if current_id in mapping_dict:
                next_ids = mapping_dict[current_id]

                # Double bonds
                if len(next_ids) == 2 and next_ids[0] == next_ids[1]:  # Check for doubles
                    doubles.add(tuple(sorted([current_id, next_ids[0]])))  # Export the double bond
                    return # Does not continue chain is double bond detected
                
                # Continue down chain
                for next_id in next_ids:
                    if next_id != prev_id: # Avoid going back to the previous node
                        stack.append((next_id, current_id))
        if ring_detected:
            rings.append(path[ring_start:] + [path[ring_start]])

        elif len(path) > 1: 
            chains.append(path)

    for key in mapping_dict.keys():
        if key not in visited:
            trace_chain(key)

    return chains, rings, list(doubles)

def getDynamicMultimerprops(rset, deltasnap, N, Rval, outputfile, outputfile2):
    print(np.shape(rset))
    r =  rset[0]
    Nsnaps = len(r.snap)
    if deltasnap>1:
        Ndatapoints = int(Nsnaps/deltasnap)+1
    else:
        Ndatapoints = int(Nsnaps)
    datawritestatus = False
    datamatrix = np.zeros((7,Ndatapoints))
    datamatrix[0][:] = np.arange(0,Nsnaps,deltasnap)
    length_df = pd.DataFrame()
    scounter = 0
    for SREF in range(0,Nsnaps,deltasnap):
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print("Snap ="+str(SREF)+" Current Time =", current_time, flush=True)

        for runidx in range(0,len(rset)):
            nsnaplocal = len(rset[runidx].snap)
            # Check if current run has the SREF'th snap. If so, collect data, if not, move on.
            if(SREF<=nsnaplocal):
                # The two variables below together represent information of each bond in the network.
                # crossboundfilid = np.array(rset[runidx].snap[SREF].crossboundfilid)
                crossboundblobid = rset[runidx].snap[SREF].crossboundblobid
                # partnerid = rset[runidx].snap[SREF].partnerid
                selfblobid = rset[runidx].snap[SREF].selfblobid
                
                datawritestatus = True
                # multimerid = rset[runidx].snap[SREF].dimerid
                partneridblob = rset[runidx].snap[SREF].partnerid
                boundid = rset[runidx].snap[SREF].boundid
                # partnermap = np.column_stack((multimerid, partneridblob))

                nHands = 3*len(rset[runidx].snap[SREF].solidcoord)
                nSol = len(rset[runidx].snap[SREF].solidcoord)

                # all_multimer_hands = [number for number in range(1,nHands+1) if number % 3 != 0]
                # monomerid = [number for number in all_multimer_hands if number not in multimerid]

                multimer_count = Counter(selfblobid)
                Multimer_1 = [num for num, freq in multimer_count.items() if freq == 1]
                Multimer_2 = [num for num, freq in multimer_count.items() if freq == 2]
                
                all_multimer_blobs = [number for number in range(1,nSol+1)]
                Multimer_0 = [number for number in all_multimer_blobs if number not in selfblobid]

                check_total = len(Multimer_0)+len(Multimer_1)+len(Multimer_2)
                if check_total != nSol:
                    print('error: M_0 + M_1 + M_2 != nSolids')

                all_actin_hands = [number for number in range(3,nHands+3) if number % 3 == 0]
                freeid = [number for number in all_actin_hands if number not in boundid]

                check_total = len(freeid)+len(crossboundblobid)
                if check_total != nSol:
                    print('error: bound + free != nSolids')

                freeid_blob = [number for number in all_multimer_blobs if number not in crossboundblobid]

                # corresponding_boundid_blob = [number/3 for number in boundid]
                bound_multimer_0 = [number for number in Multimer_0 if number in crossboundblobid]
                bound_multimer_1 = [number for number in Multimer_1 if number in crossboundblobid]
                bound_multimer_2 = [number for number in Multimer_2 if number in crossboundblobid]
                
                # corresponding_freeid_blob = [number/3 for number in freeid]
                free_multimer_0 = [number for number in Multimer_0 if number not in crossboundblobid]
                free_multimer_1 = [number for number in Multimer_1 if number not in crossboundblobid]
                free_multimer_2 = [number for number in Multimer_2 if number not in crossboundblobid]

                check_total_b = len(bound_multimer_0)+len(bound_multimer_1)+len(bound_multimer_2)
                check_total_f = len(free_multimer_0)+len(free_multimer_1)+len(free_multimer_2)
                check_total_t = check_total_b + check_total_f
                if check_total_b != len(crossboundblobid):
                    print('error: sum(bound_multimers) != crossboundblobid')
                if check_total_f != len(freeid_blob):
                    print('error: sum(free_multimers) != freeid_blob')
                if check_total_t != nSol:
                    print('error: sum(multimers) != nSolids')
                
                partnermap_blob_dict = defaultdict(list)

                for key, value in zip(selfblobid, partneridblob):
                    partnermap_blob_dict[key].append(value)

                partnermap_blob_dict = dict(partnermap_blob_dict)
                
                chains, rings, doubles = find_chains_rings_and_doubles(partnermap_blob_dict)

                chain_lengths = [len(chain) for chain in chains]
                ring_lengths = [len(set(ring)) for ring in rings]
                double_lengths = [len(double) for double in doubles]

                for double in doubles:
                    if len(double) != 2:
                        print('Warning: Double has length != 2')

                all_lengths = chain_lengths + ring_lengths + double_lengths
                length_frequencies = Counter(all_lengths)
                length_df_temp = pd.DataFrame(length_frequencies.items(), columns=['Length', f'{scounter}'], dtype='Int64')
                if length_df.empty:
                    length_df = length_df_temp
                else:
                    length_df = pd.merge(length_df, length_df_temp, on='Length', how='outer')

        if(datawritestatus):
            # Does not look across multimer bond
            datamatrix[1][scounter]= len(bound_multimer_0)
            datamatrix[2][scounter]= len(bound_multimer_1)
            datamatrix[3][scounter]= len(bound_multimer_2)
            datamatrix[4][scounter]= len(free_multimer_0)
            datamatrix[5][scounter]= len(free_multimer_1)
            datamatrix[6][scounter]= len(free_multimer_2)

            length_df.fillna(0, inplace=True)
            length_df.set_index('Length', inplace=True)
            length_df.sort_index(inplace=True)

        scounter = scounter + 1 
        
    #WRITE DATA TO FILE
    print('Saving in file named '+outputfile,flush=True)
    #Nfpc = Number of filaments per crosslinker
    pd.DataFrame(datamatrix, index=['Time', 'bound_multimer_0', 'bound_multimer_1', 'bound_multimer_2', 
                                    'free_multimer_0', 'free_multimer_1', 'free_multimer_2']).to_csv(outputfile)
    print('Saving in file named '+outputfile2,flush=True)
    pd.DataFrame(length_df).to_csv(outputfile2)

def analyzetrajectory(fpathvar,N, dirname, outfilename):
    f=4
    print(dirname,flush=True)
    meshvar=ic.icosphere(f,1.0);
    # vertices = meshvar[0]
    triangles = meshvar[1]
    Ntriangles = np.shape(triangles)[0]
    unitnormal = ic.generateunitnormals(meshvar,1.0)
    Nlists = ic.generateNeighborList(meshvar)
    print('mesh created', flush=True)
    #Generate
    deltasnap = 1
    dirlist = []
    dirlist.append(dirname)
    r = readcytosimtraj(fpathvar+dirname+'/','filonly')
    getfracoccupiedtimeseries(N, dirname, meshvar, Nlists, unitnormal, Ntriangles, r, deltasnap, outfilename)

def frontend_FO_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/FO_'+foldername+'_'+dirname+'.csv'
    print(foldername, flush=True)
    analyzetrajectory(fpathvar,N, dirname, outputfile)
    print("The End....", flush=True)

def frontend_cross_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print(args, flush=True)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/Crosslink_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)   
    rset = []
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getvaspprops(rset, 1, N, Rval,outputfile)
    print("The End....", flush=True)

def frontend_cross_LpdVASP_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print(args, flush=True)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/Crosslink_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    ratio1 = args[2]
    print(foldername, flush=True)   
    rset = []
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getvasppropsLpdVASP(rset, 1, N, Rval,outputfile, ratio1)
    print("The End....", flush=True)

def frontend_cross_DynamicDimer_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print(args, flush=True)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/Crosslink_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)   
    print(foldername, flush=True)   
    rset = []
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getDynamicDimerprops(rset, 1, N, Rval, outputfile)
    print("The End....", flush=True)

def frontend_cross_DynamicMultimer_set_Rval_repid(N, Rval, repid, *args):
    dirname = 'R_'+str(Rval)+'_r_'+str(repid)
    print(args, flush=True)
    foldername = args[0]
    fpathvar = args[1]+'/'+foldername+'/'  
    outputfile = args[1]+'outputfiles/Crosslink_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    outputfile2 = args[1]+'outputfiles/Length_Freq_'+foldername+'_R_'+str(Rval)+'_r_'+str(repid)+'.csv'
    print(foldername, flush=True)   
    rset = []
    print('Reading trajectory '+dirname,flush=True)
    r = readcytosimtraj(fpathvar+dirname+'/')
    print('Trajectory read..',flush=True)
    rset.append(r[0])
    print('Read repid '+str(repid),flush=True)
    getDynamicMultimerprops(rset, 1, N, Rval, outputfile, outputfile2)
    print("The End....", flush=True)