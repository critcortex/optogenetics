/* Created by Language version: 6.2.0 */
/* NOT VECTORIZED */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "scoplib.h"
#undef PI
 
#include "md1redef.h"
#include "section.h"
#include "md2redef.h"

#if METHOD3
extern int _method3;
#endif

#undef exp
#define exp hoc_Exp
extern double hoc_Exp();
 
#define _threadargscomma_ /**/
#define _threadargs_ /**/
 	/*SUPPRESS 761*/
	/*SUPPRESS 762*/
	/*SUPPRESS 763*/
	/*SUPPRESS 765*/
	 extern double *getarg();
 static double *_p; static Datum *_ppvar;
 
#define t nrn_threads->_t
#define dt nrn_threads->_dt
#define gkbar _p[0]
#define ik _p[1]
#define zinf _p[2]
#define gk _p[3]
#define m _p[4]
#define z _p[5]
#define cai _p[6]
#define minf _p[7]
#define mexp _p[8]
#define zexp _p[9]
#define Dm _p[10]
#define Dz _p[11]
#define _g _p[12]
#define _ion_cai	*_ppvar[0]._pval
#define _ion_ik	*_ppvar[1]._pval
#define _ion_dikdv	*_ppvar[2]._pval
 
#if MAC
#if !defined(v)
#define v _mlhv
#endif
#if !defined(h)
#define h _mlhh
#endif
#endif
 static int hoc_nrnpointerindex =  -1;
 /* external NEURON variables */
 extern double celsius;
 /* declaration of user functions */
 static int _hoc_alp();
 static int _hoc_bet();
 static int _hoc_rate();
 static int _hoc_state();
 static int _mechtype;
extern int nrn_get_mechtype();
 static _hoc_setdata() {
 Prop *_prop, *hoc_getdata_range();
 _prop = hoc_getdata_range(_mechtype);
 _p = _prop->param; _ppvar = _prop->dparam;
 ret(1.);
}
 /* connect user functions to hoc names */
 static IntFunc hoc_intfunc[] = {
 "setdata_KC3", _hoc_setdata,
 "alp_KC3", _hoc_alp,
 "bet_KC3", _hoc_bet,
 "rate_KC3", _hoc_rate,
 "state_KC3", _hoc_state,
 0, 0
};
#define alp alp_KC3
#define bet bet_KC3
 extern double alp();
 extern double bet();
 /* declare global and static user variables */
#define ek ek_KC3
 double ek = -85;
 /* some parameters have upper and lower limits */
 static HocParmLimits _hoc_parm_limits[] = {
 0,0,0
};
 static HocParmUnits _hoc_parm_units[] = {
 "ek_KC3", "mV",
 "gkbar_KC3", "mho/cm2",
 "ik_KC3", "mA/cm2",
 0,0
};
 static double delta_t = 1;
 static double m0 = 0;
 static double v = 0;
 static double z0 = 0;
 /* connect global user variables to hoc */
 static DoubScal hoc_scdoub[] = {
 "ek_KC3", &ek_KC3,
 0,0
};
 static DoubVec hoc_vdoub[] = {
 0,0,0
};
 static double _sav_indep;
 static void nrn_alloc(), nrn_init(), nrn_state();
 static void nrn_cur(), nrn_jacob();
 
static int _ode_count();
 /* connect range variables in _p that hoc is supposed to know about */
 static char *_mechanism[] = {
 "6.2.0",
"KC3",
 "gkbar_KC3",
 0,
 "ik_KC3",
 "zinf_KC3",
 "gk_KC3",
 0,
 "m_KC3",
 "z_KC3",
 0,
 0};
 static Symbol* _ca_sym;
 static Symbol* _k_sym;
 
static void nrn_alloc(_prop)
	Prop *_prop;
{
	Prop *prop_ion, *need_memb();
	double *_p; Datum *_ppvar;
 	_p = nrn_prop_data_alloc(_mechtype, 13, _prop);
 	/*initialize range parameters*/
 	gkbar = 0.08;
 	_prop->param = _p;
 	_prop->param_size = 13;
 	_ppvar = nrn_prop_datum_alloc(_mechtype, 3, _prop);
 	_prop->dparam = _ppvar;
 	/*connect ionic variables to this model*/
 prop_ion = need_memb(_ca_sym);
 nrn_promote(prop_ion, 1, 0);
 	_ppvar[0]._pval = &prop_ion->param[1]; /* cai */
 prop_ion = need_memb(_k_sym);
 	_ppvar[1]._pval = &prop_ion->param[3]; /* ik */
 	_ppvar[2]._pval = &prop_ion->param[4]; /* _ion_dikdv */
 
}
 static _initlists();
 static void _update_ion_pointer(Datum*);
 _KC3_reg() {
	int _vectorized = 0;
  _initlists();
 	ion_reg("ca", -10000.);
 	ion_reg("k", -10000.);
 	_ca_sym = hoc_lookup("ca_ion");
 	_k_sym = hoc_lookup("k_ion");
 	register_mech(_mechanism, nrn_alloc,nrn_cur, nrn_jacob, nrn_state, nrn_init, hoc_nrnpointerindex, 0);
 _mechtype = nrn_get_mechtype(_mechanism[1]);
     _nrn_thread_reg(_mechtype, 2, _update_ion_pointer);
  hoc_register_dparam_size(_mechtype, 3);
 	hoc_register_cvode(_mechtype, _ode_count, 0, 0, 0);
 	hoc_register_var(hoc_scdoub, hoc_vdoub, hoc_intfunc);
 	ivoc_help("help ?1 KC3 /home/sjarvis1/git/optogenetics/models/purkinje_miyasho/x86_64/KC3.mod\n");
 hoc_register_limits(_mechtype, _hoc_parm_limits);
 hoc_register_units(_mechtype, _hoc_parm_units);
 }
static int _reset;
static char *modelname = "BK calcium-activated potassium current";

static int error;
static int _ninits = 0;
static int _match_recurse=1;
static _modl_cleanup(){ _match_recurse=1;}
static rate();
static state();
 
static int  state (  )  {
   rate ( _threadargscomma_ v , cai ) ;
   m = m + mexp * ( minf - m ) ;
   z = z + zexp * ( zinf - z ) ;
   
/*VERBATIM*/
	return 0;
  return 0; }
 
static int _hoc_state() {
  double _r;
   _r = 1.;
 state (  ) ;
 ret(_r);
}
 
double alp (  _lv , _lca )  
	double _lv , _lca ;
 {
   double _lalp;
 _lalp = 400.0 / ( _lca * 1000.0 ) ;
   
return _lalp;
 }
 
static int _hoc_alp() {
  double _r;
   _r =  alp (  *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
double bet (  _lv )  
	double _lv ;
 {
   double _lbet;
 _lbet = 0.11 / exp ( ( _lv - 35.0 ) / 14.9 ) ;
   
return _lbet;
 }
 
static int _hoc_bet() {
  double _r;
   _r =  bet (  *getarg(1) ) ;
 ret(_r);
}
 
static int  rate (  _lv , _lca )  
	double _lv , _lca ;
 {
   double _la , _lb ;
 _la = alp ( _threadargscomma_ _lv , _lca ) ;
   zinf = 1.0 / ( 1.0 + _la ) ;
   zexp = ( 1.0 - exp ( - dt / 10.0 ) ) ;
   _lb = bet ( _threadargscomma_ _lv ) ;
   minf = 7.5 / ( 7.5 + _lb ) ;
   mexp = ( 1.0 - exp ( - dt * ( 7.5 + _lb ) ) ) ;
    return 0; }
 
static int _hoc_rate() {
  double _r;
   _r = 1.;
 rate (  *getarg(1) , *getarg(2) ) ;
 ret(_r);
}
 
static int _ode_count(_type)int _type; { hoc_execerror("KC3", "cannot be used with CVODE");}
 extern void nrn_update_ion_pointer(Symbol*, Datum*, int, int);
 static void _update_ion_pointer(Datum* _ppvar) {
   nrn_update_ion_pointer(_ca_sym, _ppvar, 0, 1);
   nrn_update_ion_pointer(_k_sym, _ppvar, 1, 3);
   nrn_update_ion_pointer(_k_sym, _ppvar, 2, 4);
 }

static void initmodel() {
  int _i; double _save;_ninits++;
 _save = t;
 t = 0.0;
{
  m = m0;
  z = z0;
 {
   rate ( _threadargscomma_ v , cai ) ;
   m = minf ;
   z = zinf ;
   }
  _sav_indep = t; t = _save;

}
}

static void nrn_init(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 v = _v;
  cai = _ion_cai;
 initmodel();
 }}

static double _nrn_current(double _v){double _current=0.;v=_v;{ {
   ik = gkbar * m * z * z * ( v - ek ) ;
   }
 _current += ik;

} return _current;
}

static void nrn_cur(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; double _rhs, _v; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
  cai = _ion_cai;
 _g = _nrn_current(_v + .001);
 	{ double _dik;
  _dik = ik;
 _rhs = _nrn_current(_v);
  _ion_dikdv += (_dik - ik)/.001 ;
 	}
 _g = (_g - _rhs)/.001;
  _ion_ik += ik ;
#if CACHEVEC
  if (use_cachevec) {
	VEC_RHS(_ni[_iml]) -= _rhs;
  }else
#endif
  {
	NODERHS(_nd) -= _rhs;
  }
 
}}

static void nrn_jacob(_NrnThread* _nt, _Memb_list* _ml, int _type){
Node *_nd; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml];
#if CACHEVEC
  if (use_cachevec) {
	VEC_D(_ni[_iml]) += _g;
  }else
#endif
  {
     _nd = _ml->_nodelist[_iml];
	NODED(_nd) += _g;
  }
 
}}

static void nrn_state(_NrnThread* _nt, _Memb_list* _ml, int _type){
 double _break, _save;
Node *_nd; double _v; int* _ni; int _iml, _cntml;
#if CACHEVEC
    _ni = _ml->_nodeindices;
#endif
_cntml = _ml->_nodecount;
for (_iml = 0; _iml < _cntml; ++_iml) {
 _p = _ml->_data[_iml]; _ppvar = _ml->_pdata[_iml];
 _nd = _ml->_nodelist[_iml];
#if CACHEVEC
  if (use_cachevec) {
    _v = VEC_V(_ni[_iml]);
  }else
#endif
  {
    _nd = _ml->_nodelist[_iml];
    _v = NODEV(_nd);
  }
 _break = t + .5*dt; _save = t;
 v=_v;
{
  cai = _ion_cai;
 { {
 for (; t < _break; t += dt) {
 error =  state();
 if(error){fprintf(stderr,"at line 52 in file KC3.mod:\n:	gk = gkbar*m*z*z\n"); nrn_complain(_p); abort_run(error);}
 
}}
 t = _save;
 } }}

}

static terminal(){}

static _initlists() {
 int _i; static int _first = 1;
  if (!_first) return;
_first = 0;
}
