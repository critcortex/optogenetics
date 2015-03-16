#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," CaE.mod");
    fprintf(stderr," CaEcvode.mod");
    fprintf(stderr," CaP2.mod");
    fprintf(stderr," CaP2cvode.mod");
    fprintf(stderr," CaT.mod");
    fprintf(stderr," CaTcvode.mod");
    fprintf(stderr," CalciumP.mod");
    fprintf(stderr," K23.mod");
    fprintf(stderr," K23cvode.mod");
    fprintf(stderr," KA.mod");
    fprintf(stderr," KAcvode.mod");
    fprintf(stderr," KC3.mod");
    fprintf(stderr," KC3cvode.mod");
    fprintf(stderr," KD.mod");
    fprintf(stderr," KDcvode.mod");
    fprintf(stderr," KM.mod");
    fprintf(stderr," KMcvode.mod");
    fprintf(stderr," Kh.mod");
    fprintf(stderr," Khcvode.mod");
    fprintf(stderr," Khh.mod");
    fprintf(stderr," Khhcvode.mod");
    fprintf(stderr," Leak.mod");
    fprintf(stderr," NaF.mod");
    fprintf(stderr," NaFcvode.mod");
    fprintf(stderr," NaP.mod");
    fprintf(stderr," NaPcvode.mod");
    fprintf(stderr, "\n");
  }
  _CaE_reg();
  _CaEcvode_reg();
  _CaP2_reg();
  _CaP2cvode_reg();
  _CaT_reg();
  _CaTcvode_reg();
  _CalciumP_reg();
  _K23_reg();
  _K23cvode_reg();
  _KA_reg();
  _KAcvode_reg();
  _KC3_reg();
  _KC3cvode_reg();
  _KD_reg();
  _KDcvode_reg();
  _KM_reg();
  _KMcvode_reg();
  _Kh_reg();
  _Khcvode_reg();
  _Khh_reg();
  _Khhcvode_reg();
  _Leak_reg();
  _NaF_reg();
  _NaFcvode_reg();
  _NaP_reg();
  _NaPcvode_reg();
}
