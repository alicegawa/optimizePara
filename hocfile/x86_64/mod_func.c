#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _pGPeA_k_reg(void);
extern void _preNoisyI_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," pGPeA_k.mod");
    fprintf(stderr," preNoisyI.mod");
    fprintf(stderr, "\n");
  }
  _pGPeA_k_reg();
  _preNoisyI_reg();
}
