import numpy as np
import math
from time import time
import datetime

start = time()
today = datetime.date.today()
Ry2eV = 13.6056981
hb = 1.054e-34
fh = 6.625e-34
q = 1.602e-19
mo = 9.11e-31
t = hb*hb/2.0/mo/q
a0 = [0]*14
str_fac = [[0]*14]*6
mat = ["Si", "Ge", "Sn", "GaP", "GaAs", "AlSb", "InP", "GaSb", "InAs", "InSb", "ZnS", "ZnSe", "ZnTe", "CdTe"]

Num_semi = 14
NN_neighbours = 10
N_pw = 137
N_bands = 12
N_ek = 40
ig2_mag = [0, 3, 4, 8, 11, 12, 16, 19, 20, 24];
Lwork = 2*N_pw - 1;
print("--- Date = %s --- Runtime = %s seconds ---" % (today,(time() - start)))
