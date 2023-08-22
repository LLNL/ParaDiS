#pragma once

#ifndef _PDS_GPU_ERROR_H
#define _PDS_GPU_ERROR_H

#include "cuda_portability.h"

class GPU_Error_Vector
{
   public :
     __cuda_hdev__  GPU_Error_Vector(void) {}
     __cuda_hdev__ ~GPU_Error_Vector()     {}

     __cuda_host__ void Allocate   (const int nerr);

     __cuda_host__ void Reset      (void);
     __cuda_host__ void Get_Errors (void);

     __cuda_hdev__ void Set_Error  (const int tndx, const int code);
};

#endif  // _PDS_GPU_ERROR_H
