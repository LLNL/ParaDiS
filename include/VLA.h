#pragma once

#ifndef _PDS_VLA_H
#define _PDS_VLA_H

/*--------------------------------------------------------------------------
 *	V3.h Declarations for 1x3 vector class
 *------------------------------------------------------------------------*/

#define VLA_VECTOR_BYTES(vtype,var,vlen)       ((vlen)*sizeof(vtype))
#define VLA_MATRIX_BYTES(mtype,var,rows,cols)  (((rows)*sizeof(mtype *)) + ((rows)*(cols)*sizeof(mtype)))

#define VLA_VECTOR(vtype,var,vlen,ptr)         vtype *var = (vtype *) ptr; ptr+=(vlen)*sizeof(vtype);
#define VLA_VECTOR_ZERO(vtype,var,vlen)        memset(var,0,(vlen)*sizeof(vtype));

#define VLA_MATRIX(mtype,var,rows,cols,ptr)    mtype **var = (mtype **) ptr; ptr+=(rows)*sizeof(mtype *); \
                                               for (int i=0; (i<(rows)); i++) \
                                               { var[i]=(mtype *) ptr; ptr+=(cols)*sizeof(mtype); }

#define VLA_MATRIX_ZERO(mtype,var,rows,cols)   { unsigned char *tmp = ((unsigned char *) var) + (rows)*sizeof(mtype *); \
                                                 memset(tmp,0,(rows)*(cols)*sizeof(mtype));                              }

#endif  // _PDS_VLA_H
