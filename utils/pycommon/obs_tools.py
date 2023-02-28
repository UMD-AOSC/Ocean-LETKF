import numpy as np
import os
from scipy.io import FortranFile
import os

def get_index(nrec):
    if nrec == 9:
        ix={"elem": 0, "rlon":1, "rlat":2, "rlev":3, "odat":4, "oerr":5, "ohx":6, "oqc":7, "obhr":8}
    else:
        raise RunTimeError("nrec={}. Unrecognized format. Exit...".format(nrec))
        sys.exit(112)
    return ix

def read_obs(fname,dtype='<f4'):
    fname = os.path.abspath(fname)
    f = FortranFile(fname,'r')
    nobs = 0
    while True:
        try:
            buf = f.read_reals(dtype=dtype)
            nobs+=1
            if nobs ==1:
                nrec = buf.size
                ix = get_index(nrec)
        except:
            print("totally {} obs from file: {}".format((nobs,nrec), fname))
            break
    f.close()

    buf2d = np.zeros((nobs,nrec),dtype=float)
    f = FortranFile(fname,'r')
    for i in range(nobs):
        buf2d[i,:] = f.read_reals(dtype=dtype)
    f.close()

    return buf2d, ix

if __name__ == '__main__':
    obs, ix = read_obs("obs01002.dat")
    print(obs[0,:])
    print(obs[-1,:])
