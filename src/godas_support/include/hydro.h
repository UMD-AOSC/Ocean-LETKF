#ifndef _HYDRO_H_
#define _HYDRO_H_

#define G0 980.0

extern float atg( float, float, float );

extern float density( float, float, float );

extern double dbl_density( float, float, float );

extern void geopan( float*, float*, float*, float*, float*, int, int );

extern float grav( float );

extern float press( float, float );

extern float theta( float, float, float, float );

extern float tempa( float, float, float, float );

extern float zeta( float, float );

#endif
