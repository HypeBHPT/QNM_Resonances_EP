#include "1+1_Free_Boundaries.h"

// -------------------------------------------------------------------------------
void allocate_derivs(derivs *V, int nl, int nh)
{	
	(*V).X   = dvector(nl, nh);
	(*V).X1  = dvector(nl, nh);
	(*V).X11 = dvector(nl, nh);
	(*V).Y   = dvector(nl, nh);
	(*V).Y1  = dvector(nl, nh);
	(*V).Y11 = dvector(nl, nh);
}
// -------------------------------------------------------------------------------
void allocate_derivs_NULL(derivs *V)
{	
	(*V).X   = NULL;
	(*V).X1  = NULL;
	(*V).X11 = NULL;
	(*V).Y   = NULL;
	(*V).Y1  = NULL;
	(*V).Y11 = NULL;
}
// -------------------------------------------------------------------------------
void fill0_derivs(derivs V, int nl, int nh)
{	
	fill0_dvector(V.X,   nl, nh);
	fill0_dvector(V.X1,  nl, nh);
	fill0_dvector(V.X11, nl, nh);
	fill0_dvector(V.Y,   nl, nh);
	fill0_dvector(V.Y1,  nl, nh);
	fill0_dvector(V.Y11, nl, nh);
}
// -------------------------------------------------------------------------------
void copy_derivs(derivs Vout, derivs Vin, int nl, int nh)
{	
	copy_dvector(Vout.X,   Vin.X,   nl, nh);
	copy_dvector(Vout.X1,  Vin.X1,  nl, nh);
	copy_dvector(Vout.X11, Vin.X11, nl, nh);
	copy_dvector(Vout.Y,   Vin.Y,   nl, nh);
	copy_dvector(Vout.Y1,  Vin.Y1,  nl, nh);
	copy_dvector(Vout.Y11, Vin.Y11, nl, nh);
}
// -------------------------------------------------------------------------------
void free_derivs(derivs *V, int nl, int nh)
{	
	free_dvector((*V).X,   nl, nh);
	free_dvector((*V).X1,  nl, nh);
	free_dvector((*V).X11, nl, nh);
	free_dvector((*V).Y,   nl, nh);
	free_dvector((*V).Y1,  nl, nh);
	free_dvector((*V).Y11, nl, nh);
}
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
// -------------------------------------------------------------------------------
void allocate_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch)
{	
	(*V).X   = dmatrix(nrl, nrh, ncl, nch);
	(*V).X1  = dmatrix(nrl, nrh, ncl, nch);
	(*V).X11 = dmatrix(nrl, nrh, ncl, nch);
	(*V).Y   = dmatrix(nrl, nrh, ncl, nch);
	(*V).Y1  = dmatrix(nrl, nrh, ncl, nch);
	(*V).Y11 = dmatrix(nrl, nrh, ncl, nch);
}
// -------------------------------------------------------------------------------
void fill0_cheb_derivs(cheb_derivs V, int nrl, int nrh, int ncl, int nch)
{	
	fill0_dmatrix(V.X,   nrl, nrh, ncl, nch);
	fill0_dmatrix(V.X1,  nrl, nrh, ncl, nch);
	fill0_dmatrix(V.X11, nrl, nrh, ncl, nch);
	fill0_dmatrix(V.Y,   nrl, nrh, ncl, nch);
	fill0_dmatrix(V.Y1,  nrl, nrh, ncl, nch);
	fill0_dmatrix(V.Y11, nrl, nrh, ncl, nch);
}
// -------------------------------------------------------------------------------
void copy_cheb_derivs(cheb_derivs Vout, cheb_derivs Vin, int nrl, int nrh, int ncl, int nch)
{	
	copy_dmatrix(Vout.X,   Vin.X,   nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.X1,  Vin.X1,  nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.X11, Vin.X11, nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.Y,   Vin.Y,   nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.Y1,  Vin.Y1,  nrl, nrh, ncl, nch);
	copy_dmatrix(Vout.Y11, Vin.Y11, nrl, nrh, ncl, nch);
}
// -------------------------------------------------------------------------------
void free_cheb_derivs(cheb_derivs *V, int nrl, int nrh, int ncl, int nch)
{	
	free_dmatrix((*V).X,   nrl, nrh, ncl, nch);
	free_dmatrix((*V).X1,  nrl, nrh, ncl, nch);
	free_dmatrix((*V).X11, nrl, nrh, ncl, nch);
	free_dmatrix((*V).Y,   nrl, nrh, ncl, nch);
	free_dmatrix((*V).Y1,  nrl, nrh, ncl, nch);
	free_dmatrix((*V).Y11, nrl, nrh, ncl, nch);
}
// -------------------------------------------------------------------------------
