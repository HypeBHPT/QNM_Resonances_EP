#include "1+1_Free_Boundaries.h"

void Get_Differentiation_Matrices(parameters *par)
{
	int iDom, N0 = (*par).N0;
	
	for(iDom=0; iDom<nDoms; iDom++){
		int N1 = (*par).N1[iDom];
		(*par).Dx1[iDom] = dmatrix(0, N1, 0, N1);
		Chebyshev_Differentiation_Matrix(N1, (*par).Dx1[iDom], (*par).grid_x1[iDom]);
	}
	(*par).Dx0 = dmatrix(0, N0, 0, N0);
	Chebyshev_Differentiation_Matrix(N0, (*par).Dx0, (*par).grid_x0);
}
// -------------------------------------------------------------------------------
void Free_Differentiation_Matrices(parameters *par)
{
	int iDom, N0 = (*par).N0;
	
	for(iDom=0; iDom<nDoms; iDom++){
		int N1 = (*par).N1[iDom];
		free_dmatrix((*par).Dx1[iDom], 0, N1, 0, N1);
	}
	free_dmatrix((*par).Dx0, 0, N0, 0, N0);
}

