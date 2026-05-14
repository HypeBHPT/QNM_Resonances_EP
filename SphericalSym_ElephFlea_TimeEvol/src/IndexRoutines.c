#include "1+1_Free_Boundaries.h"

// -------------------------------------------------------------------------------
int Index_i(parameters par, int iDom, int i)
{
	int N1 = par.N1[iDom], i1 = i;
	
	if (i1 < 0)  i1 = -i1;
	if (i1 > N1) i1 = 2*N1-i1;
	
	return i1;
}
// -------------------------------------------------------------------------------
int Index_j(parameters par, int iDom, int j)
{
	int N0 = par.N0, j1 = j;
	
	if (j1 < 0)  j1 = -j1;
	if (j1 > N0) j1 = 2*N0-j1;
	
	return j1;
}
// -------------------------------------------------------------------------------
int Index_slice(parameters par, int iDom, int iField, int i1)
{	
	int mA = par.mA[iDom], ii = Index_i(par, iDom, i1);
	return  iField + nFields*(mA + ii);
}
// -------------------------------------------------------------------------------
int Index_field1D_slice(parameters par, int iField)
{	
	int indx = par.ntotal_2Dslice + iField;
	return  indx;
}
// -------------------------------------------------------------------------------
int Index_field1D(parameters par, int iField, int i0){
	int indx, ntotal_2D=par.ntotal_2D, N0=par.N0;

	indx = ntotal_2D + iField + nFields1D*(N0-i0);

	return indx;
}
// -------------------------------------------------------------------------------
int Index(parameters par, int iDom, int iField, int i1, int i0)
{	// mA[iDom] = Sum[ nA[jDom] , {jDom, 0, iDom-1} ]; mA[0] = 0

	int mA = par.mA[iDom], N0 = par.N0,
		ii = Index_i(par, iDom, i1), jj = Index_j(par, iDom, i0);

	return  iField + nFields*(mA + ii + par.mA[nDoms]*(N0-jj));
}
// -------------------------------------------------------------------------------
void Get_Index_Arrays(parameters *par)
{	// mA[iDom] = Sum[ nA[jDom] , {jDom, 0, iDom-1} ]; mA[0] = 0
	// mA[iDom] gives total number of points up domain iDom-1
	// e.g. 2 Domains: N0=4 (n0=5), N1=3 (n1=4)
	// mA[0] = 0 (No points as there are no prvious domains)
	// mA[1] = 5 (Number of points in domain 0)
	// mA[2] = 9 (Number of points up to domain 1 --- Note there are two points at same location, i.e. the domain boundary)
	
	int iDom, N0=(*par).N0, n0=N0+1, N1;

	(*par).Nmax = maximum2((*par).N0, (*par).N1[0]);
	
	(*par).mA[0] = 0;	
	for(iDom=1; iDom<=nDoms; iDom++){
		(*par).Nmax = maximum2((*par).Nmax, (*par).N1[iDom]);
		(*par).mA[iDom] = (*par).mA[iDom-1]+((*par).N1[iDom-1]+1);
	}

	(*par).ntotal_2Dslice   = nFields*(*par).mA[nDoms]; //Total number of points at a time constant slice for 1+1D fields
	(*par).ntotal_slice   = (*par).ntotal_2Dslice + nFields1D;  //Total number of points at a time constant slice for 1+1D fields and 1D fields

	(*par).ntotal_2D     = (*par).ntotal_slice*n0; // Total number of points for 1+1D fields
	(*par).ntotal_1D   = nFields1D*n0; //Total number of points for 1D fields in time direction
	(*par).ntotal = n0*(*par).ntotal_slice + (*par).ntotal_1D; // Total number of points for 1+1D fields and 1D fields
	(*par).Ntotal = (*par).ntotal - 1;


	(*par).grid_points_x0  = dvector(0, N0);
	(*par).angle_incr_x0 = get_grid_angle_incremente((*par).grid_x0, N0);
	
	int i;
	for(i=0; i<=N0; i++){
		double x = get_grid_point((*par).grid_x0, i, N0);
		(*par).grid_points_x0[i]=x;		
	}

	for(iDom=0; iDom<nDoms; iDom++){
		N1 = (*par).N1[iDom];
		(*par).grid_points_x1[iDom]  = dvector(0, N1);
		(*par).angle_incr_x1[iDom] = get_grid_angle_incremente((*par).grid_x1[iDom], N1);

		for(i=0; i<=N1; i++){
			double x = get_grid_point((*par).grid_x1[iDom], i, N1);
			(*par).grid_points_x1[iDom][i]=x;		
		}
	}
	
	return;
}
// -------------------------------------------------------------------------------
void Get_Indices_From_Indx(parameters par, int indx, int *iDom, int *iField, int *i1, int *i0)
{
	int aux, N0=par.N0, N0mj;
	
	N0mj = indx/(nFields*par.mA[nDoms]);
	*i0   = N0-N0mj;
	
	if( *i0<=N0 ){
		aux      = indx - N0mj*nFields*par.mA[nDoms];
		*iDom    = nDoms-1;
		while(aux < nFields*par.mA[*iDom]) 
			*iDom -= 1;
		aux     -= nFields*par.mA[*iDom];
		*i1       = aux/nFields;
		*iField  = aux - (*i1)*nFields;
	}

	return;
}
// -------------------------------------------------------------------------------
void Get_Indices_From_Indx_slice(parameters par, int indx_slice, int *iDom, int *iField, int *i1)
{
	int aux;
	
		*iDom    = nDoms-1;
		while(indx_slice < nFields*par.mA[*iDom]) 
			*iDom -= 1;
		aux      = indx_slice - nFields*par.mA[*iDom];
		*i1       = aux/nFields;
		*iField  = aux - (*i1)*nFields;

		return;
}
// -------------------------------------------------------------------------------
void Free_Index_Arrays(parameters *par)
{	
	int N0=(*par).N0, N1, iDom;
	free_dvector((*par).grid_points_x0, 0, N0);

	for(iDom=0; iDom<nDoms; iDom++){
		N1 = (*par).N1[iDom];
		free_dvector((*par).grid_points_x1[iDom], 0, N1);
	}
	
	return;
}
// -------------------------------------------------------------------------------


