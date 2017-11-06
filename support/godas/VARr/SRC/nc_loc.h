#ifndef _NC_LOC_H_
#define _NC_LOC_H_

typedef struct {                             /* dimension */
               char name[MAX_NC_NAME];
               size_t size;
               } NCDIM;

typedef struct {                             /* variable */
               char name[MAX_NC_NAME];
               nc_type type;
               int ndims;
               int dims[MAX_VAR_DIMS];
               int natts;
               } NCVAR;

typedef struct {                             /* attribute */
               int var;
               char name[MAX_NC_NAME];
               nc_type type;
               int len;
               void *val;
               } NCATT;
#endif
