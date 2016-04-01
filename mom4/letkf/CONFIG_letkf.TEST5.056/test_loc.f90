PROGRAM test_loc




USE vars_letkf,  ONLY: var_local
USE letkf_local, ONLY: obs_local

INTEGER,      :: ij,ilev,mlev,nobstotal
REAL(r_size), :: var_local(nid_obs)
REAL(r_size), :: hdxf(nobstotal,nbv)
REAL(r_size), :: rdiag(nobstotal)
REAL(r_size), :: rloc(nobstotal)
REAL(r_size), :: dep(nobstotal)
INTEGER,      :: nobsl


!STEVE: need to change localization for obs below mlev - added to input arguments (3/25/2016)
ij = 1
ilev = 1
mlev = 10
CALL obs_local(ij,ilev,mlev,var_local(n,:),hdxf,rdiag,rloc,dep,nobsl,nobstotal)


END PROGRAM test_loc
