#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," exp2i.mod");
    fprintf(stderr," h.mod");
    fprintf(stderr," kadist.mod");
    fprintf(stderr," kaprox.mod");
    fprintf(stderr," kdrca1.mod");
    fprintf(stderr," na3s.mod");
    fprintf(stderr," naxn.mod");
    fprintf(stderr," netstims.mod");
    fprintf(stderr," nmdanet.mod");
    fprintf(stderr, "\n");
  }
  _exp2i_reg();
  _h_reg();
  _kadist_reg();
  _kaprox_reg();
  _kdrca1_reg();
  _na3s_reg();
  _naxn_reg();
  _netstims_reg();
  _nmdanet_reg();
}
