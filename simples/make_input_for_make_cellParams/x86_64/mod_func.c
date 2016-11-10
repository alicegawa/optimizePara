#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;

extern void _hh_k_reg(void);
extern void _pGPeA_mine_reg(void);

void modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," hh_k.mod");
    fprintf(stderr," pGPeA_mine.mod");
    fprintf(stderr, "\n");
  }
  _hh_k_reg();
  _pGPeA_mine_reg();
}
