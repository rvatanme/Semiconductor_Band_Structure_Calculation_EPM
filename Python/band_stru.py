import numpy as np
import math
from time import time
import datetime
import os
import cmath
import scipy.linalg as la

def generate_Hamil(Inp, dk1, Ig_Vec0, pseu, ek, H_EPM):
    pi = 4.0*math.atan(1.0)
    k0 = 2.0*pi/Inp['ao']
    N_pw = 137
    
    for i in range(N_pw):
        Hrow = []
        gx_i = Ig_Vec0[i]['igx_vec']
        gy_i = Ig_Vec0[i]['igy_vec']
        gz_i = Ig_Vec0[i]['igz_vec']
        for j in range(N_pw):
            gx_j = Ig_Vec0[j]['igx_vec']
            gy_j = Ig_Vec0[j]['igy_vec']
            gz_j = Ig_Vec0[j]['igz_vec']
            if (i == j):
                gkx2 = (gx_i+dk1['dkx'])*(gx_i+dk1['dkx'])
                gky2=(gy_i+dk1['dky'])*(gy_i+dk1['dky'])
                gkz2=(gz_i+dk1['dkz'])*(gz_i+dk1['dkz'])
                Hreal= Inp['t']*k0*k0*(gkx2+gky2+gkz2)
                Hrow.append(complex(Hreal, 0.0))
                #if (i < 5):
                #    H1.append(complex(Hreal, 0.0))
                #    print("i=%s\tt=%s\tk0=%s\tgkx2=%s\tgky2=%s\tgkz2=%s"%(i,Inp['t'],k0,gkx2,gky2,gkz2))
                #H_EPM[i][j].imag=0.0
            else:
                dgx=gx_i-gx_j
                dgy=gy_i-gy_j
                dgz=gz_i-gz_j
                idgx=round(dgx)
                idgy=round(dgy)
                idgz=round(dgz)
                id2=idgx*idgx+idgy*idgy+idgz*idgz
                tau_fac=2.00*pi/8.00
                if (id2==3):
                    V_s=pseu['Vs3']
                    V_a=pseu['Va3']
                elif (id2==4):
                    V_s=pseu['Vs4']
                    V_a=pseu['Va4']
                elif (id2==8):
                    V_s=pseu['Vs8']
                    V_a=pseu['Va8']
                elif (id2==11):
                    V_s=pseu['Vs11']
                    V_a=pseu['Va11']
                else:
                    V_s=0.00
                    V_a=0.00
                ct=math.cos((dgx+dgy+dgz)*tau_fac)
                st=math.sin((dgx+dgy+dgz)*tau_fac)
                Hreal= V_s*ct
                Himag= V_a*st
                #if (i == 0 and j < 5):
                #    print("(i,j)=(%s, %s)\tdgx=%s\tdgy=%s\tdgz=%s\ttau_fac=%s\tHreal=%s"%(i,j,dgx,dgy,dgz,tau_fac,Hreal))
                Hrow.append(complex(Hreal, Himag))
                #H_EPM[i][j] = complex(Hreal, Himag);
                #H_EPM[i][j].imag= Himag;
        H_EPM.append(Hrow)
    
    '''for i in range(5):
        for j in range(5):
            print ("%s\t"%H_EPM[i][j].real)
        print ("\n")'''
    #for i in range(5):
    #    print(H1[i].real)
    eigvals, eigvecs = la.eig(H_EPM)
    Real = eigvals.real
    Real.sort()
    for i in range(8):
        ek.append(Real[i])
    #print (ek)
    return (0)
    

start = time()
today = datetime.date.today()
Ry2eV = 13.6056981
hb = 1.054e-34
fh = 6.625e-34
q = 1.602e-19
mo = 9.11e-31
t = hb*hb/2.0/mo/q
a0 = []
str_fac = []
mat = ["Si", "Ge", "Sn", "GaP", "GaAs", "AlSb", "InP", "GaSb", "InAs", "InSb", "ZnS", "ZnSe", "ZnTe", "CdTe"]

Num_semi = 14
NN_neighbours = 10
N_pw = 137
N_bands = 12
N_ek = 40
ig2_mag = [0, 3, 4, 8, 11, 12, 16, 19, 20, 24];
Lwork = 2*N_pw - 1;

f1 = open("a0.txt", "r")
lines=f1.readlines()
for i in range(len(lines)):
    a0.append(float(lines[i]))  
f1.close()

f1 = open("str_fac.txt", "r")
lines = f1.readlines()
for i in range(14):
    srow = []
    for j in range(6):
        srow.append(float(lines[i].split()[j]))
    str_fac.append(srow)
f1.close()

pseu = {'Vs3' : 0.0, 'Vs4' : 0.0,  'Vs8': 0.0, 'Vs11': 0.0, 'Va3' : 0.0, 'Va4' : 0.0, 'Va8' : 0.0, 'Va11' : 0.0}
dk1 = {'dkx': 0.0, 'dky' : 0.0, 'dkz' : 0.0}
Ig_vec = {'igx_vec': 0, 'igy_vec' : 0, 'igz_vec' : 0}
Ig_Vec0 = []
Inp = {'ao' : 0.0, 't' : 0.0}
i = 0

icount = 0
for j in range(NN_neighbours):
    ig2 = ig2_mag[j];
    a = math.sqrt(ig2)
    ia = math.floor(a)
    ic = 0
    for igx in range(-ia, ia +1):
        igx2 =  igx*igx
        for jgy in range(-ia, ia + 1):
            jgy2 = jgy*jgy
            igxy2 = igx2 + jgy2
            for kgz in range(-ia, ia + 1):
                kgz2 = kgz*kgz
                igxyz2 = igxy2 + kgz2
                if (igxyz2 == ig2):
                    ds = {}                    
                    ds['igx_vec'] = igx
                    ds['igy_vec'] = jgy
                    ds['igz_vec'] = kgz
                    Ig_Vec0.append(ds)
                    ic += 1
                    icount += 1
    #print ("%s %s %s %s\n" % (ig2,ia,ic,icount))


#for i in range(137):
#    print ("i=%s\tigx=%s\tigy=%s\tigz=%s"%(i,Ig_Vec0[i]['igx_vec'],Ig_Vec0[i]['igy_vec'],Ig_Vec0[i]['igz_vec']))


i = 0

while (i < 1):
    filename = mat[i] + "_py.txt"
    #print (filename)
    f1 = open(filename, "w")
    f1.write("Band structure of %s\n"%mat[i])
    f1.close()
    
    pseu['Vs3']=str_fac[i][0]*Ry2eV;
    pseu['Vs4']=0.0;
    pseu['Vs8']=str_fac[i][1]*Ry2eV;
    pseu['Vs11']=str_fac[i][2]*Ry2eV;

    pseu['Va3']=str_fac[i][3]*Ry2eV;
    pseu['Va4']=str_fac[i][4]*Ry2eV;
    pseu['Va8']=0.0;
    pseu['Va11']=str_fac[i][5]*Ry2eV;

    Inp['ao'] = a0[i]*1.0e-10
    Inp['t'] = t

    # Symmetry path L->Gamma
    ik_curr=0
    sumx=0.0
    X0=0.5
    Y0=0.5
    Z0=0.5
    
    for ik in range(1, N_ek+1):
        H_EPM = []
        ek = []
        dk1['dkx']=0.500*(1.00 - (ik-1)/N_ek)
        dk1['dky']=0.500*(1.00 - (ik-1)/N_ek)
        dk1['dkz']=0.500*(1.00 - (ik-1)/N_ek)
        generate_Hamil(Inp, dk1, Ig_Vec0, pseu, ek, H_EPM)
        ekm.append(ek)
    
        
    # Dumping data for Symmetry path L->Gamma
    myfile = open(filename, "a")
    for lh in range(1,N_ek+1):
        myfile.write(str(ik_curr + lh) + "\t")
        for kh in range(8):
            myfile.write(str(ekm[lh][kh]) + "\t")
        myfile.write("\n")
    myfile.close()
    
    
    i += 1

print("--- Date = %s --- Runtime = %s seconds ---" % (today,(time() - start)))
