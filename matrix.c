/* matrix.c: contains functions for manipulating matrices 
             and calculation RUS frequencies 
*/

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>
#include <complex.h>
#include "matrix.h"


// General mathematical functions ///////////
static double sa;
#define SQR(a) ((sa = (a)) == 0 ? 0 : sa * sa)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double delta (int a, int b){double tmp = (a == b);return tmp;}
double pythag (double a, double b)
	{
	double aa, bb; aa = fabs (a); bb = fabs (b);
	if (aa > bb) return (aa * sqrt (1.0 + SQR(bb/aa)));
	else 
		{
		if (bb == 0.0) return 0; else return (bb * sqrt (1.0 + SQR(aa/bb)));		
		}
	}

// Comparison functions for qsort
int compare_doubles (const void *a, const void *b)
	{
	const double *da = (const double *) a;
	const double *db = (const double *) b;
	return (*da > *db) - (*da < *db);
	}

int compare_indfrth (const void *f1, const void *f2)
	{
    const index_frth *ff1 = (const index_frth *)f1;
    const index_frth *ff2 = (const index_frth *)f2;
    return ( ((* ff1).fr > (* ff2).fr) - ((* ff1).fr < (* ff2).fr) );                 
	}

int compare_solution (const void *f1, const void *f2) 
	{
    const solution *ff1 = (const solution *)f1;
    const solution *ff2 = (const solution *)f2;
    return ( ((* ff1).chi2 > (* ff2).chi2) - ((* ff1).chi2 < (* ff2).chi2) );                 
	}
///////////////////////////////////////////////////////


////////////////////////////////////////////
indices getind (int index, int nmax) // Decompose an index into i,n,m,l.
	{
	indices a;
	a.i = index%4;
	index /= 4;
	a.n = index % ++nmax;
	index /= nmax;
	a.m = index % nmax;
	a.l = index / nmax;
	return a;
	}

// Core integral over volume for RUS calculation
double integ (int xexp, int yexp, int zexp, dimensions a, double *dfact)  // a contains log of dimensions, and geometry.
	{
		int prod;
		int xes,yes,zes,xp,yp,zp;
		double integral,p,r,v,t,mm,xs,ys,zs,odidr;
		prod = (xexp+1) * (yexp+1) * (zexp+1); // Used to determine parity.
		if (prod %2 == 0 && a.geom!=4) 
		{    
			return 0;  // if any of the exponents are odd, integral is zero (except potato geometry)
		}
		else 
		{
			switch (a.geom)
			{
				case 1: //RPR geometry
					integral=(exp((xexp+1)*a.logdx+(yexp+1)*a.logdy+(zexp+1)*a.logdz)*8./prod); // RPR geometry
					break;
					
				case 0: //Cylinder geometry
					p=dfact[xexp-1 +2];   // Notice the +2 needed to address arguments ranging from -1 to n
					r=dfact[yexp-1 +2];
					t=dfact[xexp+yexp+2  +2];
					integral=4*PI*(exp((xexp+1)*a.logdx+(yexp+1)*a.logdy+(zexp+1)*a.logdz))*p*r / ((zexp+1)*t);
					break;
					
				case 4: //Potato geometry
					p=dfact[xexp-1 +2];
					r=dfact[yexp-1 +2];
					v=dfact[zexp-1 +2];
					t=dfact[xexp+yexp+zexp+3 +2];
					
					//Sum on quadrants
					integral=0.;
					if((xexp%2)+(yexp%2)+(zexp%2)>1)
						{
							mm=1.;
						}
					else
						{
						mm=PI/2.;
						}
					
					for (xp=-1;xp<=1;xp=xp+2)
					{
						xs=(a.dx*(xp+1)+a.dxm*(1-xp))/2.; // x+ or x- depending on quadrant
						xes=(1-(1-xp)*(xexp%2));        // sign for quadrant contribution // -1 only if xp=-1 and xexp is odd
						for (yp=-1;yp<=1;yp=yp+2)
						{
							ys=(a.dy*(yp+1)+a.dym*(1-yp))/2.;
							yes=(1-(1-yp)*(yexp%2)); 
							for (zp=-1;zp<=1;zp=zp+2)
							{
							zs=(a.dz*(zp+1)+a.dzm*(1-zp))/2.;
							zes=(1-(1-zp)*(zexp%2)); 
							integral += xes*yes*zes*(pow(xs,xexp+1) * pow(ys,yexp+1) * pow(zs, zexp+1));
							}
						}
					}
					integral *= mm*p*r*v/t;
					break;
					
				case 5: //Sphere geometry
					p=dfact[xexp-1 +2];
					r=dfact[yexp-1 +2];
					v=dfact[zexp-1 +2];
					t=dfact[xexp+yexp+zexp+3 +2];
					integral=4*PI*exp((xexp+1)*a.logdx+(yexp+1)*a.logdy+(zexp+1)*a.logdz)*p*r*v/t;
					break;
					
				case 6: //Octahedron geometry
					p=dfact[100+xexp];
					r=dfact[100+yexp];
					v=dfact[100+zexp];
					t=dfact[100+xexp+yexp+zexp+3];
					integral=exp((xexp+1)*a.logdx+(yexp+1)*a.logdy+(zexp+1)*a.logdz)*8.*p*r*v/t;
					break;    
					
				case 8: //Hollow cylinder geometry. Makes difference between two cylinders.

					p=dfact[xexp-1 +2];   // Notice the +2 needed to address arguments ranging from -1 to n
					r=dfact[yexp-1 +2];
					t=dfact[xexp+yexp+2  +2];
					odidr=a.OD/a.ID;
					//   integral=((pow(a.OD,xexp+yexp+2)- pow(a.ID,xexp+yexp+2)) * pow(a.dz, zexp+1))*PI*(p*r*4./(zexp+1)/t);
					integral=( ((pow(odidr,xexp+yexp+2)-1.))*pow(a.ID,xexp+yexp+2) * pow(a.dz, zexp+1))*PI*(p*r*4./(zexp+1)/t);
					// The formula above works. Alternatives lead to non-positive defined matrix for EK.
					break;
					
				ruserror("Error: shape is not yet defined. Choices are 0=Cyl, 1=RPR, 4=Potato, 5=Sphere, 6=Octahedron, and 8=Hollow Cyl\n");
				
			}
			return integral;
		}
	}

// Printing and debugging functions ///////////////////////////////////////////////////////////////////
/* Print out the error message and exit. */
void ruserror (char *error_text) {fprintf(stderr, "matrix: %s\n", error_text);exit (1);}
/* Print blank line */
void println(void) {printf("\n");}
/* Print a matrix */
void print_matrix (double **mat, int dim)
	{
	int i, j; /* Iterators */
	for (i = 1; i <= dim; i++) { for (j = 1; j <= dim; j++) printf ("%g\t ", mat[i][j]); printf("\n"); }
	}
/* Print a vector with elements of type double */
void print_vector (double *vec, int dim)
	{int i; for (i = 1; i <= dim; i++){printf ("%g\n", vec[i]);}}
void print_vector_h (double *vec, int dim)
	{int i; for (i = 1; i <= dim; i++){printf ("%8.6g ", vec[i]);}printf("\n");}
/* Print a vector with elements of type integers */
void print_ivector (int *vec, int dim)
	{ int i; for (i = 1; i <= dim; i++) { printf ("%i\n", vec[i]);} }
/* Print the no to no+n elements of a vector */
void print_output (double *vec, int no, int n)
	{ int i; for (i = no; i <= no+n; i++) printf("%8.6f\n", vec[i]); }
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// General matrix and table operations //////////////////////////////////////////////////////////
/* Return a pointer to a two dimensional square array of 'dim' rows
and columns. The indices of rows and column start from 1 */
double ** matrix (int dim)
	{
	int rw = dim;
	int cl = dim;
	double **mat;
	int i; /* loop iterator */

	/* allocate an array of pointers to rows*/
	mat = (double **) malloc ((rw + 1) * sizeof(int *));
	if (! mat) 
		{
		printf("Dimensions %i %i\n",rw,cl);
		ruserror ("Failed to allocate memory in matrix()");
		}

	/* allocate rows */
	mat[1] = (double *) malloc ((rw * cl + 1) * sizeof(double));
	if (! mat[1]) ruserror ("Failed to allocate memory in matrix()");
	memset (mat[1], 0, ((rw * cl + 1) * sizeof (double))); /* sanity for Valgrind */

	for (i = 2; i <= rw; i++) mat[i] = mat[i-1] + cl;

	return mat;
	}  

/* Free the pointer to the two dimensional array 'mat' */
void free_matrix(double **mat)
{
 free (mat[1]); 
 free (mat);
}

double ** table (int dim1,int dim2)
	{
	int rw = dim1;
	int cl = dim2;
	double **mat;
	int i; /* loop iterator */

	/* allocate an array of pointers to rows*/
	mat = (double **) malloc ((rw + 1) * sizeof(int *));
	if (! mat) ruserror ("failed to allocate memory in matrix()");

	/* allocate rows */
	mat[1] = (double *) malloc ((rw * cl + 1) * sizeof(double));
	if (! mat[1]) ruserror ("failed to allocate memory in matrix()");
	memset (mat[1], 0, ((rw * cl + 1) * sizeof (double))); /* sanity for Valgrind */

	for (i = 2; i <= rw; i++) mat[i] = mat[i-1] + cl;

	return mat;
	}  

void free_table(double **mat)
	{
	free (mat[1]); 
	free (mat);
	}

/* Return a vector of dimension 'dim' with indices ranging from 1 to dim */
double *vector (int dim)
	{
	double *vec = (double *) malloc ((dim + 1) * sizeof (double));
	memset (vec, 0, ((dim + 1) * sizeof (double))); /* sanity for Valgrind */
	return vec;
	}

/* Free the pointer to the vector 'vec' */
void free_vector (double *vec) {free (vec);}
	
/* Return a freq vector of dimension 'dim' with indices ranging from 1 to dim */
index_frth *indfrth (int dim)
	{
	index_frth *object = (index_frth *) malloc ((dim + 1) * sizeof (index_frth));
	memset (object, 0, ((dim+1) * sizeof(double)));
	//memset (vec, 0, ((dim + 1) * sizeof (double))); /* sanity for Valgrind */
	return object;
	}

void free_indfrth(index_frth *indfrth) {free (indfrth); }

/* Return a solution vector of dimension 'dim' with indices ranging from 1 to dim */
solution *sol (int dim)
	{
	solution *object = (solution *) malloc ((dim + 1) * sizeof (solution));
	memset (object, 0, ((dim+1) * sizeof(double)));
	//memset (vec, 0, ((dim + 1) * sizeof (double))); /* sanity for Valgrind */
	return object;
	}

void free_solution(solution *sol) {free (sol); }

/* Return a int vector of dimension 'dim' with indices ranging from 1 to dim */
int *ivector (int dim)
	{
	int *vec = (int *) malloc ((dim +1) * sizeof (int));
	memset (vec, 0, ((dim + 1) * sizeof (int)));
	return vec;
	}

void free_ivector (int *vec) {free (vec);}

/* Calculate the length of a vector */
double lengthgrad(double *vec, double *vec2,int dim)
	{
	int i;
	double l=0;
	for (i = 1; i <= dim; i++) l += SQR(vec[i]-vec2[i]);
	return sqrt(l);
	}


/* Multiply a matrix 'mat' to the right by a vector 'vec' and return
the resulting vector */
double *mul_matrix_vec (double **mat, double *vec, int dim)
	{  
	double *mvec = vector (dim);
	int i, j; 
	for (i = 1; i <= dim; i++)
	{
		for (j = 1; j <= dim; j++) mvec[i] += (mat[i][j] * vec[j]);
	}
	return mvec;
	}

/* Multiply a matrix 'mat' by a scalar 'sca' and return the resulting matrix */
void mul_matrix_sca (double **mat, double sca, int dim)
	{  
	int i, j; 
	for (i = 1; i <= dim; i++)
	{
	for (j = 1; j <= dim; j++) mat[i][j] *= sca;
	}
	}

/* Returns the transposed of a square matrix of dimension 'dim' */
double **transpose_matrix (double **mat, int dim)
	{
	int i, j;
	double **tmat = matrix(dim);
	for (i = 1; i <= dim; i++) for (j = 1; j <= dim; j++) tmat[i][j] = mat[j][i];
	return tmat;
	}

/* Returns the product of the matrices 'mata' and 'matb' */
double **mul_matrix_matrix (double **mata,double **matb,int dim)
	{
	int i, j, k;
	double **matp = matrix (dim);
	for (i = 1; i <= dim; i++) for (j = 1; j <= dim; j++) for (k = 1; k <= dim; k++) matp[i][j] += mata[i][k] * matb[k][j];
	return matp;
	}

	
// Matrix operations for eigenvalue search.	
/* Tridiagonalizes the matrix 'a'. */
void householder (double **a, 
int n,       /* dimension of a */
double *d,   /* returns diagonal elements */
double *e)   /* returns off-diagonal elements */
	{
	int i, j, k, l;  
	double scale, hh, h, g, f;

	for (i = n; i >= 2; i--)
		{
		l = i - 1;
		h = 0.0;  /* = |x.x| */
		scale = 0.0;

		if (l > 1) 
			{
			for (k = 1; k <= l; k++)
			scale += fabs(a[i][k]);
			if (scale == 0.0)
				e[i] = a[i][l];
			else 
				{
				for (k = 1; k <= l; k++)  /* scale a's */
					{
					a[i][k] /= scale;
					h += a[i][k] * a[i][k];
					}

				f = a[i][l];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale * g;
				h -= f*g;
				a[i][l] = f - g;
				f = 0.0;

				for (j = 1; j <= l; j++)
					{
					g = 0.0;
					for (k = 1; k <=j; k++)
					g += a[j][k] * a[i][k];
					for (k = j+1; k <= l; k++)
					g += a[k][j] * a[i][k];
					e[j] = g/h;
					f += e[j] * a[i][j];
					}
				
				hh = f /(h + h);
				
				for (j = 1; j <= l; j++)
					{
					f = a[i][j];
					e[j] = g = e[j] - hh*f;
					
					for (k = 1; k <= j; k++)
					a[j][k] -= (f*e[k] + g*a[i][k]);
					}
				}
			}
		else
			{
			e[i] = a[i][l];
			}
		d[i] = h;
		}

	e[1] = 0.0;
	for (i = 1; i <= n; i++) d[i] = a[i][i];
	}

/* Returns eigenvalues in d */
void tri_eig (double *d, /* input diagonal elements */
	double *e, /* input off-diagonal elements */
	int n)    /* input dimension */
	{
	int m, l, iter, i;
	double s, r, p, g, f, dd, c, b;

	for (i = 2; i <= n; i++) e[i-1] = e[i];  /* shift the off-diagonal elements up */

	e[n] = 0.0;

	for (l = 1; l <= n; l++)
	{
	iter = 0;
	do 
	{
	for (m = l; m <= n-1; m++)
	{
		dd = fabs (d[m]) + fabs (d[m+1]);
		if ((double)(fabs(e[m]) + dd) == dd ) 
		break;
	}
	if (m != 1)
	{

		if (iter++ == 1000) {
		ruserror("Too many iterations in tri_eig");}
		g = (d[l+1] - d[l]) / (2.0 * e[l]);
		r = pythag (g, 1.0);
		g = d[m] - d[l] + e[l] / (g + SIGN(r,g));
		s = c = 1.0;
		p = 0.0;
		
		for (i=m-1; i >= l; i--)
		{
		f = s * e[i];
		b = c * e[i];
		e[i+1] = (r = pythag(f, g));
		if (r == 0.0)
		{
		d[i+1] -= p;
		e[m]  = 0.0;
		break;
		}
		s = f / r;
		c = g / r;
		g = d[i+1] - p;
		r = (d[i] - g) * s + 2.0 * c * b;
		d[i+1] = g + (p = s*r);
		g = c * r - b;
	}
	
	if (r == 0.0 && i >= l) continue;
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
	}
	} while (m != l);
	}
	}

/* Returns Cholesky decomposition. */
double **cholesky (double **mat, int dim)
	{
	double **a = matrix(dim); 
	int i, j, k;
	double sum;

	for (i = 1; i <= dim; i++)
	{
	for (j = i; j <= dim; j++)
	{
	for (sum = mat[i][j], k = i-1; k >= 1; k--)
	{
		sum -= a[i][k] * a[j][k];
	}

	if (i == j) /* Diagonal elements */
	{
		if (sum <= 0.0)
		{
		printf ("i,j: %i sum: %g\n", i,sum);
		//sum = -sum;
		ruserror ("The input matrix is not positive definite.");
		}
		a[i][i] = sqrt(sum);
		//printf("a[%d][%d] = %f\n", i, i, a[i][i]);
		}
	else /* Off-diagonal elements */
	{
		a[j][i] = sum / a[i][i];
		//printf("a[%d][%d] = %f\n", j, i, a[j][i]);
	}
	}
	}  
	return a;
	}

/* Make the 'C' matrix for Cholesky decomposition*/
void make_c (double **a, double **l, int dim)
	{
	double sum;
	int i, j, k;
	double **lt = transpose_matrix(l, dim);

	for (i = 1; i <= dim; i++)
	{
	for (j = 1; j <= dim; j++)
	{
	sum = a[i][j];
	for (k = j-1; k >= 1; k--)
	{
		sum -= a[i][k] * lt[k][j];
	}
	a[i][j] = sum / lt[j][j];
	}
	}

	for (j = 1; j <= dim; j++)
	{
	for (i = 1; i <= dim; i++)
	{
	sum = a[i][j];
	for (k = i-1; k >= 1; k--)
	{
		sum -= l[i][k] * a[k][j];
	}
	a[i][j] = sum / l[i][i];
	}
	}
	free_matrix(lt);
	}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/*
All parameters are read and stored in two vectors.
param for double parameters
[1]-[21] elastic constants. upon read same order as for rus.exe. **************************

[30]-[32] dimensions, in cm

[33]-[35] Euler Angles

[36] relative change in chisq
[37] length of gradient vector
[38] flambda
[39] density, computer generated, in mg/cm3
[40] mass, in g
[41] RMS
[53] Tracking parameter

[60] Evolutionary pressure for gen-alg. Number >1. Default is 3.
[60] For fitting: will contain chi2 stoppng value.
[61] Mixing parameters for gen-alg. - allows extrapolation. Number 0<x<0.5. Default is 0.2.
[62] Mutation rate. p/nc (number of constants). Default is p=1.

iparm for integer parameters **********************************************************
[1]-[21] 0 if parameter is free; 1 if it has to be fixed.
[30]-[32] computer generated - dimensions are free or fixed
[33]-[35] constraints on angles

[36] number of moduli to fit, sets the space group. 
	2 = isotropic
	3 = cubic
	5 = hexagonal
	6 = tetragonal [406 or 407]
	7 = trigonal   [306 or 307]
	9 = orthorombic
	13 = monoclinic
	21 = triclinic
[37] numdim 0 = fixed dimensions, >0 = free dimensions
[38] icalc  0 = fit, !=0 = number of lines to calculate or choose mode.
[39] ideriv 0 = no deritavies for simylation, 1 = output derivatives
[40] iteration counter in curfit
[41] nterms = number of fit parameters
[42] nmax, order of the polynomials
[43] dim, computer generated, dimension of the full matrix
[44] n, number of experimental peaks that are being read.
[45] shape. 5 = sphere; 1 = RPR; 0 = cylinder; 8 = hollow cylinder. Give OD, ID, height
[46] has mirror planes. Number between 0 and 3.

[60] For fitting: 10^-p = chi2 stopping criteration. Default is -5. 
[60] Genetic population size for gen-alpg. Default it 500. 
[61] Genetic algorithm: Number of generations for gen alg. Default it 10. 
*/
//////////////////////////////////////////////////////////////////////////////////////////////////////////

void read_input (char *filename, comments *COMM, int *iparm, double *param,long *fex,long *w)
	{
	FILE *fp;
	char *line = NULL;
	char buffer[256],null[256];
	char *token;
	int nlines;
	int i,j,nargs_line1;
	bool cont;
	double *ft = vector(1000);
	double *wt = vector(1000);

	if ((fp = fopen (filename, "r")) == NULL) ruserror ("cannot open input file");
	for(i=0;i<32;i++) {(*COMM).block[i]=7;}

	fgets(buffer,256,fp);  /* read line 1 - Should be title card*/
	if (buffer[0]!='#')  {sprintf((*COMM).comment[0],"# %s",buffer);}
	else                 {strcpy((*COMM).comment[0],buffer);} 

	token =strtok(buffer,":");
	while(token!=NULL) {param[53]=atof(token);token =strtok(NULL,":");}
	printf("Tracking parameter: %g\n",param[53]);
	
	//printf("Comment mark in %s %s %i\n",buffer,(*COMM).comment[0],(buffer[0]=='#') );
	i=0; j=0;
	cont=true;


	while(cont)	
	{
		fgets(buffer,256,fp);  
		if (buffer[0]=='#')
		{ 
			i++; strcpy((*COMM).comment[i],buffer); (*COMM).block[i]=1;                     //This is the first block of lines.
		}
		/* read line 2 - Shape[45] - Bravais class[36] - Max poly order[42]  Fit or calc nlines[38] */
		/* Mirror [46] , Mass [40], RMS [41] (not used), Calculate derivatives (0/1). */
	else		
	{
		nargs_line1=sscanf (buffer, "%i %i %i %i %i %lg %lg %i %i %i %lg %lg %lg",
		&iparm[45],&iparm[36],&iparm[42],&iparm[38],&iparm[46],&param[40],&param[41],&iparm[39],
		&iparm[60],&iparm[61],&param[60],&param[61],&param[62]);
		cont=false;
	}
	}
	
	// Additions for supplemental control parameters.

	if(iparm[38]==0)
		{
		if(nargs_line1>8)            // Sets chi2 stopping parameter for fitting.
			{param[60]=1./pow(10.,iparm[60]);}
		else {param[60]=0.00001;}  
		if(nargs_line1<9)            // Default to rigorous error calculation.
			{iparm[61]=1.;}  
		else 
			{
			if(iparm[61]>0.){iparm[61]=1.;}else{iparm[61]=0.;} // Deactivate rigorous error calculation only if <=0 is given.
			}
		}

	if(iparm[38]==-1||iparm[38]==-5) // Simple MC or GenAlg
		{
			if(nargs_line1>8) {if(iparm[60]<100){iparm[60]=500;}} // Set default size for simple MC simulation; if MC pop too small, set 500.
			else {iparm[60]=500;}
		}
	
	if(iparm[38]==-5) // Gen Alg
			{
			// Sets controls for genetic algorithm only if all parameters are given. Otherwise default.
			if(nargs_line1<13) {iparm[61]=10;param[60]=3.;param[61]=0.2;param[62]=1.;} // If insufficient parameters, set default.
			else{																					// If sufficient parameters, catch out-of-bound parameters.
				if(iparm[61]<2){iparm[61]=10;}	// If too small number of generations is given, default to 10.
				if(param[60]<1.5){param[60]=3.;}	// If too small evoluationary pressure is given default to 3.
				if(param[61]<0.||param[61]>0.5){param[61]=0.2;}	// If too large mixing extrapolation is given, default to 0.2. 0 is acceptable.
				if(param[62]<0.||param[62]>5.){param[62]=1.;}	// If out of bound mutation rate is given, default to 1./nc. Otherwise param[62]/nc.
				}
		}
	
	if(iparm[46]<0||iparm[46]>3)	 ruserror("Error: only 0, 1, 2, or 3 mirror planes are allowed.");
	if(iparm[36]%100>10&&iparm[36]%100!=13&&iparm[36]%100!=21) ruserror("Error: allowable choices of # of c_xx are 2, 3, 5, 306-307, 406-407, 13, or 21.");

	iparm[43] = ((iparm[42]+1) * (iparm[42]+2) * (iparm[42]+3)) / 2;// Dimension of big matrix.

	cont=true;
	while(cont)	
	{
		fgets(buffer,256,fp);                                                            // Reading line with elastic constants
		if (buffer[0]=='#')
		{
			i++; strcpy((*COMM).comment[i],buffer); (*COMM).block[i]=2;                     //This is the second block of lines.
		}
	/* read line 3 - Elastic constants */
		else
		{
			param[1]=atof(strtok(buffer," "));
			/* For trigonal and tetragonal systems, Bravais codes are 306,307, 406,407, but only 6 or 7 moduli are given */
			for (j=2;j<iparm[36]%100+1;j++) param[j]=atof(strtok(NULL," "));
			for (j=iparm[36]%100+1;j<22;j++) {param[j]=0;   iparm[j]=1;} // Initialize all to zero and fixed.
			cont=false;
	}      
	}

	fgets(buffer,256,fp);                                                            // Reading constraints on elatic constants.
		
	if (buffer[0]=='#')  // 
	{
		ruserror("Error: please specify constraints on the elastic constants in input. Aborting.\n");  //This should not happen.
	}
	else
	{
		iparm[1]=atoi(strtok(buffer," "));
		/* For trigonal and tetragonal systems, Bravais codes are 36,37, 46,47, but only 6 or 7 moduli are given */  		
		for (j=2;j<iparm[36]%100+1;j++) {iparm[j]=atoi(strtok(NULL," "));}
		/* 0 = free, 1 = fixed, positive number up to 100 = variation in % for Monte Carlo approach*/            
	}

	cont=true;

	while(cont)
	{
		fgets(buffer,256,fp);  
		if (buffer[0]=='#')
		{
			i++; strcpy((*COMM).comment[i],buffer); (*COMM).block[i]=3;                     //This is the third block of lines.
			
		}
	/* read line 4 - Dimensions */
		else
		{
			if(iparm[45]==4) sscanf (buffer, "%lg %lg %lg %lg %lg %lg",&param[30],&param[31],&param[32],&param[50],&param[51],&param[52]); 
			else sscanf (buffer, "%lg %lg %lg",&param[30],&param[31],&param[32]);
			// Get xm,ym,zm dimensions for potato
			cont=false;
		}
	}

	fgets(buffer,256,fp);  
	if (buffer[0]=='#')
	{
		ruserror("Error: please specify constraint on dimensions. 0=free, 1=fixed. Aborting.\n");  //This should not happen.
	}
	/* read line 4 - Dimensions constraint. */
	else {sscanf (buffer, "%i",&iparm[37]);}

	cont=true;

	if(iparm[37]!=0) {iparm[30]=1;iparm[31]=1;iparm[32]=1;} // Dimensions are constrained
				else {iparm[30]=0;iparm[31]=0;iparm[32]=0;} // Dimensions are free
				
	while(cont)
	{
		fgets(buffer,256,fp);  
		if (buffer[0]=='#')
		{
			i++; strcpy((*COMM).comment[i],buffer); (*COMM).block[i]=4;                     //This is the fourth block of lines.
		}
	/* Read line 5 - Euler Angles */
		else
		{
			sscanf (buffer, "%lg %lg %lg",&param[33],&param[34],&param[35]);
			cont=false;
		}
	}

	fgets(buffer,256,fp);  
	if (buffer[0]=='#')
	{
		ruserror("Error: please specify constraint on dimensions. 0=free, 1=fixed. Aborting.\n");  //This should not happen.
	}
	/* Read line 5 - Euler Angle constraints. */
	else
	{
		sscanf (buffer, "%i %i %i",&iparm[33],&iparm[34],&iparm[35]);
	}
	//printf("Dimensions: %g %g %g %g %g %g \n",param[10],param[30],param[11],param[31],param[12],param[32]);		
	param[39] = 1000*param[40]/(param[30]*param[31]*param[32]); // Get density
	switch(iparm[45])
	{
		case 1: param[39] *= 1; break;       // RP
		case 0: param[39] *= (4/PI); break;  // Elliptical cylinder
		case 4: param[39] = 1000*param[40]/(param[30]+param[50])/(param[31]+param[51])/(param[32]+param[52])*(6/PI); break;
		case 5: param[39] *= (6/PI); break;  // Spheroid // V=4/3*pi/8=pi/6
		case 6: param[39] *= (6); break;  // Octahedron // V=1/6*d1*d2*d3
		case 8: param[39] = 1000.*param[40]/(param[32]*PI/4*(param[30]*param[30]-param[31]*param[31])); break; // Hollow cylinder. OD**2 - ID**2.
		ruserror("Error: not a supported shape parameter. Use 0=Cyl, 1=RPR, 5=Spheroid, 6=Octahedron, or 8=Hollow Cyl.\n");
	}
	printf("Density: %g g/cc\n",param[39]/1000.);
	cont=true;
	while(cont)
	{
		fgets(buffer,256,fp);  
		if (buffer[0]=='#')
		{
			i++; strcpy((*COMM).comment[i],buffer); (*COMM).block[i]=5;                     //This is the fifth block of lines.
		}
	/* read list of frequencies */
		else	{cont=false;} // Found the first data line
	}

	j=1;
	sscanf (buffer, "%lg %s %lg", &ft[j],null,&wt[j]);
	cont=fgets(buffer,256,fp);
	while (cont)
	{
		if (buffer[0]!='#')
		{
			j++;sscanf (buffer, "%lg %s %lg", &ft[j],null,&wt[j]);
			cont=fgets(buffer,256,fp);
		}
		else
		{
			cont=false;
			i++; strcpy((*COMM).comment[i],buffer); (*COMM).block[i]=6; //This is the sixth, last block of lines.
		}
	}

	if((*COMM).block[i]==6) cont=true; // We did not yet reach EOF
	while(cont)
	{
		cont=fgets(buffer,256,fp);  
		//printf("Just read buffer %s",buffer);
		if (buffer[0]=='#') {i++; strcpy((*COMM).comment[i],buffer); (*COMM).block[i]=6;}           //Sixth block continuing.
		// If any lines are not commented out, they are ignored and delete from output.
		strcpy(buffer," ");
		//printf("==\n%s : %i\n==\n",(*COMM).comment[i],(*COMM).block[i]);
	}

	nlines=j+1;
	double *fex_loc=vector(nlines);
	double *w_loc=vector(nlines);
	for(j=1;j<nlines;j++)
	{
		fex_loc[j]=ft[j]; w_loc[j]=wt[j];     
	}

	// Checking consistency...
	// If Euler !=0 or not fixed warning if cylinder, spheroid.
	int irot;
	irot=3-iparm[33]-iparm[34]-iparm[35]+(param[33]!=0)+(param[34]!=0)+(param[35]!=0);

	if(iparm[45]!=1&&irot!=0){
	printf("\n=====================================================\n WARNING: using rotations with cylinder or spheroids \n=====================================================\n");
	}

	if(iparm[46]>0&&irot!=0){
	printf("=====================================================\n WARNING: using rotations with mirror planes         \n=====================================================\n");
	}

	if(iparm[46]>0&&iparm[45]==4){
	printf("=====================================================\n WARNING: using potato with mirror planes         \n=====================================================\n");
	}

	if(iparm[30]==0&&iparm[45]==4){
	printf("=====================================================\n WARNING: fitting shape deactivated for potato         \n=====================================================\n");
	iparm[30]=1;iparm[31]=1;iparm[32]=1;iparm[37]=1;
	}

	
	iparm[44]=nlines-1;
	fex[0] = (long) fex_loc;
	w[0] = (long) w_loc;
	free_vector(ft);
	free_vector(wt);
	if (line) free (line);
	fclose (fp);
	}

// Write rusio file (only after fitting)
void write_io (char *filename, comments *COMM, int *iparm, double *param,long *fex,long *fth, long *w,long *fsym)
{
FILE *fp;

int i = 0;
int j = 0;
bool cont;
if ((fp = fopen ("rusio.dat", "w")) == NULL) ruserror ("cannot open output file");

fprintf(fp,"%s",(*COMM).comment[0]);  // Write title line
j++;

cont=true;
while(cont)
{
    if((*COMM).block[j+1]==1)
    {
        fprintf(fp,"%s",(*COMM).comment[j]); // This is not the last line of block 1; print it.
        j++;
    }
    else
    {
        fprintf(fp,"# Shape     #-cxx     n-poly     #-calc    Mirror   Mass(g)  LastRMS  CalcDeriv\n"); 
        cont=false;
        // This is the last line of block 1. Print expected format.
    }
}

while((*COMM).block[j]==1){j++;}

fprintf(fp,"  %i         %i         %i         %i         %i        %8.6f %8.6f %i\n",
iparm[45],iparm[36],iparm[42],iparm[38],iparm[46],param[40],param[41],iparm[39]);

cont=true;        
while(cont)
{
    if((*COMM).block[j+1]==2)
    {
        fprintf(fp,"%s",(*COMM).comment[j]);  // This is not the last line of block 2; print it.
        j++;
    }
    else
    {
        fprintf(fp,"# Elastic constants (2nd line: n constraints (0/1))\n"); 
        cont=false;
        // This is the last line of block 2. Print expected format.
    }
}

while((*COMM).block[j]==2){j++;}


fprintf(fp,"  %8.6f ",param[1]);
for (i=2;i<iparm[36]%100+1;i++) fprintf(fp,"%8.6f ",param[i]);
fprintf(fp,"\n");

fprintf(fp,"  %i        ",iparm[1]);
for (i=2;i<iparm[36]%100+1;i++) fprintf(fp,"%i        ",iparm[i]);
fprintf(fp,"\n");

cont=true;
while(cont)
{
    if((*COMM).block[j+1]==3)
    {
        fprintf(fp,"%s",(*COMM).comment[j]);  // This is not the last line of block 3; print it.
        j++;
    }
    else
    {
        if(iparm[45]!=4)
            {fprintf(fp,"# Dimensions        (2nd line: 1 constraint (0/1))\n");}
        else
            {fprintf(fp,"# Dimensions x+ y+ z+ x- y- z- (2nd line: 1 constraint (0/1) - fitting not supported for potato shape)\n");}                
        cont=false;
        // This is the last line of block 3. Print expected format.
    }
}

while((*COMM).block[j]==3){j++;}

fprintf(fp,"  %8.6f %8.6f %8.6f",param[30],param[31],param[32]);
if(iparm[45]==4) {fprintf(fp," %8.6f %8.6f %8.6f",param[50],param[51],param[52]);}   
fprintf(fp,"\n  %i\n",iparm[37]);

cont=true;
while(cont)
{
    if((*COMM).block[j+1]==4)
    {
        fprintf(fp,"%s",(*COMM).comment[j]);  // This is not the last line of block 4; print it.
        j++;
    }
    else
    {
        fprintf(fp,"# Euler angles      (2nd line: 3 constraints (0/1))\n");
        cont=false;
        // This is the last line of block 4. Print expected format.
    }
}

while((*COMM).block[j]==4){j++;}

fprintf(fp,"  %8.6f %8.6f %8.6f\n",param[33],param[34],param[35]);
fprintf(fp,"  %i        %i        %i\n",iparm[33],iparm[34],iparm[35]);

cont=true;
while(cont)
{
    if((*COMM).block[j+1]==5)
    {
        fprintf(fp,"%s",(*COMM).comment[j]);  // This is not the last line of block 5; print it.
        j++;
    }
    else
    {
        fprintf(fp,"# F_ex(MHz) F_th(MHz) Weight\n");
        cont=false;
    }
}

while((*COMM).block[j]==5){j++;}

for(i=1;i<=iparm[44];i++) {fprintf(fp, "  %8.6f  %8.6f  %4.2f\n",((double *)fex[0])[i],((double *)fth[0])[i],((double *)w[0])[i]);}

while((*COMM).block[j]==6)
{
fprintf(fp,"%s",(*COMM).comment[j]); j++;
}


//fprintf(fp,"\n");

fclose (fp);
	
}

// Write rusout file (only after fitting)
void write_out (char *filename, comments *COMM, int *iparm, double *dfact, double *param,long *fex,long *fth, long *w,long *fsym,double *sigmap, long *dfdc,long *lookup)
{
FILE *fp, *su;

double **c = matrix (6);
double **s = matrix (6);
double sum; //,rms2;
double KV,KR,GV,GR,KVRH,GVRH,MU,VL,VT,Va,RHO;
bool   cont;
int i,j,k = 0;
int iline;
int m[22];
int kcount[22];
//int nparam;

int *dvar=ivector(35); // Inverse pairing of parameter to index in sigmap
j=1;
for(i=1;i<iparm[36]%100+1;i++) 
{
    if( iparm[i]!=0) {} // If an elastic constants is fixed
    else {dvar[i]=j;j++;}        // If it is free pair parameter dvar[i] with j
   }   
if(iparm[30]!=1){dvar[30]=j;dvar[31]=j+1;j=j+2;} // If dimensions are free, 2 additional parameters. 12 stays constrained by 10 and 11.
if(iparm[33]!=1){dvar[33]=j;j++;} // Euler angle psi
if(iparm[34]!=1){dvar[34]=j;j++;} // Euler angle theta
if(iparm[35]!=1){dvar[35]=j;j++;} // Euler angle phi
//for(i=1;i<16;i++) {printf("dvar[%i]=%i\n",i,dvar[i]);}

if ((fp = fopen ("rusout.dat", "w")) == NULL) ruserror ("cannot open output file");

fprintf(fp,"Resonant ultrasound spectroscopy analysis. Units cm, g, GPa, MHz.\n");
fprintf(fp,"%s",(*COMM).comment[0]); // Write title card
iline=1;
cont=true;

// In rusout file only the first block (minus last line) and last block of comments is written.
while(cont)
{   
    //printf("%s %i %i\n",(*COMM).comment[iline],iline,(*COMM).block[iline]);
    if((*COMM).block[iline+1]==1) fprintf(fp,"%s",(*COMM).comment[iline]); // Write comment up to last line of first block
    if((*COMM).block[iline]==6) fprintf(fp,"%s",(*COMM).comment[iline]); // Write last block    
    iline++;
    if(iline>31) cont=false;
}

fprintf(fp,"\n");

fprintf(fp,"Free moduli are");

switch(iparm[36])
   {
    case 2:
         m[1]=11;
         m[2]=44;
	 break;
    case 3:
         m[1]=11;
         m[2]=12;
	 m[3]=44;
	 break;
    case 5:
	m[1]=33;
	m[2]=23;
	m[3]=12;
	m[4]=44;
	m[5]=66;
	    break;
    case 406:
     m[1]=11;
	 m[2]=33;
	 m[3]=23;
	 m[4]=12;
	 m[5]=44;
	 m[6]=66;
         break;
    case 407:
     m[1]=11;
	 m[2]=33;
	 m[3]=23;
	 m[4]=12;
	 m[5]=44;
	 m[6]=66;
     m[7]=16;
         break;
    case 306:   //c11,c33,c13,c12,c44,c14
         m[1]=11;
         m[2]=33;
         m[3]=23;
         m[4]=12;
         m[5]=44;
         m[6]=14;
         break;
    case 307:   //c11,c33,c13,c12,c44,c14
         m[1]=11;
         m[2]=33;
         m[3]=23;
         m[4]=12;
         m[5]=44;
         m[6]=14;
         m[7]=25;
         break;
    case 9:   
	m[1]=11;
	m[2]=22;
	m[3]=33;
	m[4]=23;
	m[5]=13;
	m[6]=12;
	m[7]=44;
	m[8]=55;
	m[9]=66;
        break;   
    case 13:   
	m[1]=11;
	m[2]=22;
	m[3]=33;
	m[4]=23;
	m[5]=13;
	m[6]=12;
	m[7]=44;
	m[8]=55;
	m[9]=66;
    m[10]=15;
	m[11]=25;
	m[12]=35;
	m[13]=46;
	    break;   
    case 21:   
	m[1]=11;
	m[2]=22;
	m[3]=33;
	m[4]=23;
	m[5]=13;
	m[6]=12;
	m[7]=44;
	m[8]=55;
	m[9]=66;
    m[10]=14;
    m[11]=15;
	m[12]=16;
	m[13]=24;
	m[14]=25;
	m[15]=26;
	m[16]=34;
	m[17]=35;
	m[18]=36;
	m[19]=45;
    m[20]=46;
    m[21]=56;
	    
        break;   
   }

for(i=1;i<iparm[36]%100+1;i++) {if(iparm[i]==0) fprintf(fp," c%i,",m[i]);}

fprintf(fp,"\nUsing %ith order polynomials. Mass = %8.4lg g; rho = %8.6lg g/cm3\n\n",iparm[42],param[40],param[39]/1000.);
if(iparm[46]==0) fprintf(fp,"NOT ");
fprintf(fp,"Using mirror planes.\n");
fprintf(fp,"\n");

fprintf(fp,"  n    fexp(MHz)    fcalc(MHz)   %%err wt     k    i     ");

 

if(iparm[39]==1)
{
    gen_dfdc(dfdc,lookup,param,iparm,dfact,fth,fsym);
    for(i=1;i<iparm[36]%100+1;i++) { fprintf(fp,"dlnf/dlnc%i ",m[i]);}
    fprintf(fp,"sum");
} 
for(i=1;i<22;i++) {kcount[i]=1;}

fprintf(fp,"\n");
for(i=1;i<=iparm[44];i++)
    {
        sum=0;
        fprintf(fp, "%3i    %8.7f    %8.7f   %5.2f %4.2f   %1i  %3i    ", i,
                                ((double *)fex[0])[i],
                                ((double *)fth[0])[i],
                                100*((((double *)fth[0])[i]/((double *)fex[0])[i])-1.),
                                ((double *)w[0])[i],((int *)fsym[0])[i],kcount[((int *)fsym[0])[i]]);
        kcount[((int *)fsym[0])[i]]++;
// Added code here for option to calculate df/d(moduli)
        if(iparm[39]==1){
        for(j=1;j<iparm[36]%100+1;j++)
        {
                sum+= ((double *)dfdc[0])[(j-1)*iparm[44]+i];
                fprintf(fp,"%8.5f    ",((double *)dfdc[0])[(j-1)*iparm[44]+i]);
        }
        fprintf(fp,"%8.5f",sum);}
        fprintf(fp,"\n");
    }
if(iparm[39]==1){free_vector((double *)dfdc[0]);}          
fprintf(fp,"\n");

gen_elatensor(param,iparm,c);
gen_elatensor(param,iparm,s);
mul_matrix_sca (c,10.,6); // This affects also the relative error, see below.
mul_matrix_sca (s,10.,6); // S will contain the compliance tensor


for (k=1;k<=6;k++)   
  { s[k][k]= -1/s[k][k];   
   for (i=1;i<=6;i++)             
        if (i!=k) s[i][k]*=s[k][k];
    for (i=1;i<=6;i++)
        if (i!=k) 
            for (j=1;j<=6;j++)
                      if (j!=k)
                          s[i][j]+=s[i][k]*s[k][j];
    for (i=1;i<=6;i++)
       if (i!=k) 
            s[k][i]*=s[k][k];
   }
   for (i=1;i<=6;i++)        
     for (j=1;j<=6;j++)
        s[i][j]=-s[i][j];




KV=(c[1][1]+c[2][2]+c[3][3]+2*(c[1][2]+c[2][3]+c[3][1]))/9.;
KR=1./(s[1][1]+s[2][2]+s[3][3]+2*(s[1][2]+s[2][3]+s[3][1]));
GV=(c[1][1]+c[2][2]+c[3][3]-(c[1][2]+c[2][3]+c[3][1])+3*(c[4][4]+c[5][5]+c[6][6]))/15.;
GR=15./(4*(s[1][1]+s[2][2]+s[3][3])-4*(s[1][2]+s[2][3]+s[3][1])+3*(s[4][4]+s[5][5]+s[6][6]));
KVRH=(KV+KR)/2.;
GVRH=(GV+GR)/2.;
MU=(3*KVRH-2*GVRH)/(6*KVRH+2*GVRH);
RHO=param[39]*1000.;
VL=sqrt((KVRH+4*GVRH/3.)/RHO)*1000.;
VT=sqrt(GVRH/RHO)*1000.;
Va=1./cbrt((1./3)*(1/pow(VL,3.)+2/pow(VT,3.)));

fprintf(fp," Bulk modulus(GPa) KV = %6.4lg | KR = %6.4lg\n",KV,KR);
fprintf(fp,"\n");
fprintf(fp,"Shear modulus(GPa) GV = %6.4lg | GR = %6.4lg\n",GV,GR);
fprintf(fp,"\n");
fprintf(fp,"     Poisson ratio mu = %6.4lg\n",MU);
fprintf(fp,"\n");
fprintf(fp,"Speed of sound: longitudinal = %6.4lg (km/s), transverse = %6.4lg (km/s), average = %6.4lg (km/s)\n",VL,VT,Va);
fprintf(fp,"\n");
/**** Adapt for monoclinic, triclinic*/
fprintf(fp,"    c11      c22      c33      c23      c13      c12      c44      c55      c66      c14      c15      c25      c35      c16      c46\n");
fprintf(fp," %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg\n",c[1][1],c[2][2],c[3][3],c[2][3],c[1][3],c[1][2],c[4][4],c[5][5],c[6][6],c[1][4],c[1][5],c[2][5],c[3][5],c[1][6],c[4][6]);
fprintf(fp,"\n");

if ((su = fopen("summary.dat", "r")))  { fclose(su); if ((su = fopen ("summary.dat", "a")) == NULL) ruserror ("cannot open output file");}
else {if ((su = fopen ("summary.dat", "w")) == NULL) ruserror ("cannot open output file");
          fprintf(su,"# Param.    rms%%     c11      c22      c33      c23      c13      c12      c44      c55      c66      c14      c15      c25      c35      c16      c46   Poisson v_long  v_tran v_aver   dx     dy     dz     phi      theta    psi \n");            
    }


fprintf(su,"%8.5lg %8.5lg ",param[53],param[41]);
fprintf(su,"%8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg %8.5lg ",c[1][1],c[2][2],c[3][3],c[2][3],c[1][3],c[1][2],c[4][4],c[5][5],c[6][6],c[1][4],c[1][5],c[2][5],c[3][5],c[1][6],c[4][6]);
fprintf(su,"%6.4lg %6.4lg  %6.4lg %6.4lg  ",MU,VL,VT,Va);
fprintf(su,"%6.4f %6.4f %6.4f ",param[30],param[31],param[32]);
fprintf(su,"%8.5f %8.5f %8.5f\n",param[33],param[34],param[35]);

 

fclose(su);

free_matrix(c);
free_matrix(s);

fprintf(fp,"  dx(cm) dy(cm) dz(cm)    ");
fprintf(fp,"Shape set to ");
switch(iparm[45]) {
case 0:
fprintf(fp,"cylinder.\n");
break;
case 1:
fprintf(fp,"rectangular parallelepiped.\n");
break;
case 4:
fprintf(fp,"potato (x+ y+ z+ x- y- z-).\n");
break;
case 5:
fprintf(fp,"spheroid.\n");
break;
case 8:
fprintf(fp,"hollow cylinder (d1/d2 are OD/ID).\n");
break;
}
fprintf(fp," %6.4f %6.4f %6.4f",param[30],param[31],param[32]);
if(iparm[25]==4) {fprintf(fp," %6.4f %6.4f %6.4f",param[50],param[51],param[52]);}
fprintf(fp,"\n");
fprintf(fp,"\n");
fprintf(fp,"  psi      theta    phi\n");
fprintf(fp," %8.5f %8.5f %8.5f\n",param[33],param[34],param[35]);
fprintf(fp,"\n");

fprintf(fp,"loop# %i   rms error = %6.4lg %%, changed by %8.6lg %%\n",iparm[40],param[41],param[36]*100.);
fprintf(fp,"\n");

//fprintf(fp,"length of gradient vector = %lg  blamb = %lg\n",param[17],param[18]);
//fprintf(fp,"\n");

fprintf(fp,"Estimated error on free parameters - double check with error calculation tool:\n");

//rms2 = param[41]*param[41];

for(i=1;i<22;i++) 
    {
     if(iparm[i]==0) fprintf(fp," c%i      ",m[i]); 
    }  
if(iparm[30]==0) fprintf(fp," d1[cm]   d2[cm]   d3[cm]  ");
for(i=33;i<36;i++) 
    {
     switch(i){
         case 33:if(iparm[i]==0) fprintf(fp,"phi[o]     ");break;
         case 34:if(iparm[i]==0) fprintf(fp,"theta[o]   ");break;
         case 35:if(iparm[i]==0) fprintf(fp,"psi[o]     ");break;     }
    }  
fprintf(fp,"\n");

for(i=1;i<22;i++) 
    {
     //if(iparm[i]==0) fprintf(fp," %6.4f %% ",sqrt(0.01*rms2)*sigmap[dvar[i]]/param[i]);
	if(iparm[i]==0) fprintf(fp," %6.4f %% ",100.*sigmap[dvar[i]]/param[i]); 
    }  
    
double d1,d2,d3,ed1,ed2,ed3,vo;
    
if(iparm[30]==0)
    {
    d1=param[30];
    d2=param[31];
    d3=param[32];
    vo=param[40]/param[39]; // Get volume from mass/density.
    ed1=100.*sigmap[dvar[30]]/param[30];
    ed2=100.*sigmap[dvar[31]]/param[31];
    ed3=(vo/d1*ed1+vo/d2*ed2)/(vo/d3);
    fprintf(fp," %6.4f %%  %6.4f %%  %6.4f %% ",ed1,ed2,ed3);
    }   
for(i=33;i<36;i++) 
    {
             if(iparm[i]==0) fprintf(fp," %6.4f deg ",sigmap[dvar[i]]);
    }  
   
    
fprintf(fp,"\n");
fprintf(fp,"\nFor publications resulting from the use of this software, please cite J. Torres et al., JASA 151, 3547 (2022)\n");

free_ivector(dvar);
fclose (fp);
}

	
/* Generates lookup tables. Either with or without mirror planes */
void gen_lookups (long *lookup, int nmax, int mirr)
	{
	int id;  /* index to track i, l, m, n. */
	int i, l, m, n,pk;
	int kx[9],ky[9],kz[9];
	int dim = (nmax+1) * (nmax+2) * (nmax+3) / 2;
	int *kdim = ivector (9);
	int *looki = ivector (dim); /* lambda; pg 39 */
	int *lookpk = ivector (dim);  /* parity lookup table */

	for(i=1;i<9;i++) 
	{
	kdim[i]=0;
	kx[i]=i;
	ky[i]=8+((i+1)%2)*2-i;
	kz[i]=13-8*(i>4)-ky[i]; 
	// This does some magic that generates the recipe for block diagonalization.
	// kx={1,2,3,4,5,6,7,8}
	// ky={7,8,5,6,3,4,1,2}
	// kz={6,5,8,7,2,1,4,3}
	// There is a (.+4)%8 difference with Migliori
	}

	id = 0;
	for (l = 0; l < nmax+1; l++)
	{
	for (m = 0; m < nmax + 1 -l; m++)
		{
		for (n = 0; n < nmax + 1 - l -m; n++)
			{
			for (i = 1; i < 4; i++)
				{
				id++;
				looki[id] = i + 4 * (n + (nmax + 1) * (m + (nmax +1) * l));
	
				if(mirr==3)
				{ // Eight possible blocks if all mirror planes
				//pk=5-4*(l%2)+2*(m%2)+(n%2); // Parity of the triplet (p.45)
				pk=1+4*(l%2)+2*(m%2)+(n%2); // Parity of the triplet (p.45)
				// Now depending on whether i=1,2,3 (x,y,z)
				// Assign the index to one of the 8 block diagonal matrices
				if(i==1) pk=kx[pk];
				if(i==2) pk=ky[pk];
				if(i==3) pk=kz[pk];
				kdim[pk]++;
				lookpk[id] = pk; // For each id we store to which of the 8 submatrices it belongs
				}
				
				if(mirr==2)
				{ // Four possible blocks if yz and xy are mirror planes
				//pk=5-4*(l%2)+(n%2); // Parity of the triplet (p.45) -- Switch m off. 4 blocks
				pk=1+2*(l%2)+(n%2); // Parity of the triplet (p.45)
				// Now depending on whether i=1,2,3 (x,y,z)
				// Assign the index to one of the 4 block diagonal matrices
				ky[1]=3;ky[2]=4;ky[3]=1;ky[4]=2;
				kz[1]=4;kz[2]=3;kz[3]=2;kz[4]=1;
				//ky[5]=1;ky[6]=2;ky[1]=5;ky[2]=6;       
				// 3,4,7,8 are empty.
				if(i==1) pk=kx[pk];
				if(i==2) pk=ky[pk];
				if(i==3) pk=kz[pk];
				kdim[pk]++; lookpk[id] = pk; // For each id we store to which of the 8 submatrices it belongs
				}

				if(mirr==1)
				{ // Two possible blocks only if xy is a mirror planes
				pk=(n%2)+1;
				if(i==3) pk=2-(n%2);
				kdim[pk]++; lookpk[id] = pk; // For each id we store to which of the 8 submatrices it belongs
				}
				
				if(mirr==0) {kdim[1]++;lookpk[id] = 1;} // Only one block if no mirror planes
														// All polynomials are stored in one submatrix
				}
				
			}
		}
	}
   
	lookup[0]= (long) looki;
	lookup[1]= (long) lookpk;
	lookup[2]= (long) kdim;
	}

/* Generates kinetic and potential energy matrices */
void gen_ke_pe_mat (long *RA,long *lookup,int nmax,double rho,dimensions dimsample,double **c,int bdi, angles euler, double *dfact)
{

int i, j,k,l,ip,jp,kp,lp; 
int ct[4][4]; //pseudo tensor, index 1 to 4
int irot;
double d[4][4][4][4]; // rotated tensor
double rot[4][4]; // rotation matrix
double psi,theta,phi; // Euler angle - Kocks convention
double dt; //buffer for elastic tensor
int dim = (nmax+1) * (nmax+2) * (nmax+3) / 2; /* dimension of EK and EP */
double modtemp; 
dimensions logdim;

if(euler.a==0 && euler.b ==0 && euler.c ==0)
{
    irot=0;
}
else // Apply rotation. Kock's angles. Wikipedia>Euler angles>Z1Y2Z3 with 3 -> pi-3 (phi = pi -phi).
{
    irot=1;
    psi=euler.a*PI/180;
    theta=euler.b*PI/180;
    phi=euler.c*PI/180;
    //rot[row][column]
    rot[1][1]=-sin(psi)*sin(phi)-cos(psi)*cos(phi)*cos(theta);
    rot[1][2]= cos(psi)*sin(phi)-sin(psi)*cos(phi)*cos(theta);
    rot[1][3]= cos(phi)*sin(theta);

    rot[2][1]= sin(psi)*cos(phi)-cos(psi)*sin(phi)*cos(theta);
    rot[2][2]=-cos(psi)*cos(phi)-sin(psi)*sin(phi)*cos(theta);
    rot[2][3]= sin(phi)*sin(theta);

    rot[3][1]= cos(psi)*sin(theta);
    rot[3][2]= sin(psi)*sin(theta);
    rot[3][3]= cos(theta);
}


/* Generate Voigt-convention pseudo-tensor */
ct[1][1] = 1;
ct[2][2] = 2;
ct[3][3] = 3;
ct[2][3] = 4;
ct[3][2] = 4;
ct[3][1] = 5;
ct[1][3] = 5;
ct[1][2] = 6;
ct[2][1] = 6;
double **EK = matrix (((int *)lookup[2])[bdi]); /* Kinetic energy matrix */
double **EP = matrix (((int *)lookup[2])[bdi]); /* Potential energy matrix */
int iprim, jprim; //
int bdii,bdiiprim;
int factor; //
int p, q, r; //
indices indi, indiprim; //
double buff; //


if(irot==1) 
{
    
for (ip =1;ip<4;ip++) { for (jp=1;jp<4;jp++) { for (kp=1;kp<4;kp++) { for (lp=1;lp<4;lp++) {
dt=0;
for (i=1;i<4;i++) { for (j=1;j<4;j++) { for (k=1;k<4;k++){ for (l=1;l<4;l++) {

dt += rot[ip][i]*rot[jp][j]*rot[kp][k]*rot[lp][l]* c[ ct[i][j] ][ ct[k][l] ]; 

}}}}
d[ip][jp][kp][lp]=dt;
}}}}

    
}

/* Generate kinetic and potential energy matrices */
// adjust dimensions to match formula for integral
//dimsample.dx/=2.;
//dimsample.dy/=2.;
//dimsample.dz/=2.;
//Calculate log only once and feed for integral
if(dimsample.geom!=4)
    {
    logdim.dx=(dimsample.dx/2.);
    logdim.dy=(dimsample.dy/2.);
    logdim.dz=(dimsample.dz/2.);
    logdim.logdx=log(dimsample.dx/2.);
    logdim.logdy=log(dimsample.dy/2.);
    logdim.logdz=log(dimsample.dz/2.);
    }
else
    {
    logdim.dx=(dimsample.dx);
    logdim.dy=(dimsample.dy);
    logdim.dz=(dimsample.dz);
    logdim.dxm=(dimsample.dxm);
    logdim.dym=(dimsample.dym);
    logdim.dzm=(dimsample.dzm);
    }

logdim.geom=dimsample.geom;
if(logdim.geom==8){logdim.ID=(dimsample.dy/2.);logdim.OD=(dimsample.dx/2.);}

//for(bdi=1;bdi<9;bdi++)
//   {
//   printf("%g %g %g %g %g %g %g %g\n",c[1][1],c[4][4],logdim.dx,logdim.dxm,logdim.dy,logdim.dym,logdim.dz,logdim.dzm);
   bdii=0;
   bdiiprim=0;
   for (i = 1; i <= dim; i++)
      {
      indi = getind (((int *)lookup[0])[i], nmax);
      if(bdi==((int *)lookup[1])[i])
         {
         bdiiprim=0;
         bdii++;
         for (iprim = 1; iprim <= dim; iprim++)
            {
            indiprim = getind (((int *)lookup[0])[iprim], nmax);
            if(bdi==((int *)lookup[1])[iprim])
               {
               bdiiprim++;
               p = indi.l + indiprim.l;
               q = indi.m + indiprim.m;
               r = indi.n + indiprim.n;
               if(indi.i==indiprim.i)
               {
               EK[bdii][bdiiprim] = rho* integ (p, q, r, logdim, dfact);
               //printf("%i %i %.32e\n",bdii,bdiiprim,EK[bdii][bdiiprim]);
               }
               else
               {    
               EK[bdii][bdiiprim] =0;
               }    
               buff = 0;
               for (j = 1; j < 4; j++)
                  {
                  for (jprim = 1; jprim < 4; jprim++)
                     {
                     p = indi.l + indiprim.l - (j == 1) - (jprim == 1);
                     q = indi.m + indiprim.m - (j == 2) - (jprim == 2);
                     r = indi.n + indiprim.n - (j == 3) - (jprim == 3);
                     factor = ((indi.l * (j == 1) + indi.m * (j == 2) + indi.n * (j == 3)) 
                              * (indiprim.l * (jprim == 1) + indiprim.m 
                              * (jprim == 2) + indiprim.n * (jprim == 3)));
                     if(irot==0) modtemp=c[ct[indi.i][j]][ct[indiprim.i][jprim]];
                     if(irot==1) modtemp=d[indi.i][j][indiprim.i][jprim];

                     if(modtemp!=0&&factor!=0)
                        {
                        buff += (modtemp *factor * integ (p, q, r, logdim, dfact) );
                        }
                     }
                  }
                  EP[bdii][bdiiprim] = buff;
               }
            }
         }
     }

   RA[0] = (long) EK;
   RA[1] = (long) EP;
   
}

void gen_elatensor(double *param, int *iparm, double **c)
{
double *mod=vector(22);
int i;

for(i=1;i<22;i++) {mod[i]=0;}
switch(iparm[36])
   {
    case 2:   // Isotropic c11 c44
        mod[1]=param[1];
        mod[2]=mod[1];
        mod[3]=mod[1];
        mod[4]=mod[1]-2*param[2];
        mod[5]=mod[4];
        mod[6]=mod[4];
        mod[7]=param[2];
        mod[8]=mod[7];
        mod[9]=mod[7];
        break;
    case 3:   // Cubic c11 c12 c44
        mod[1]=param[1];
        mod[2]=mod[1];
        mod[3]=mod[1];
        mod[4]=param[2];
        mod[5]=mod[4];
        mod[6]=mod[4];
        mod[7]=param[3];
        mod[8]=mod[7];
        mod[9]=mod[7];
        break;
    case 5:   // Hexagonal c33 c23 c12 c44 c66
        mod[3]=param[1];
        mod[4]=param[2];
        mod[6]=param[3];
        mod[7]=param[4];
        mod[9]=param[5];
        mod[1]=2*mod[9]+mod[6];
        mod[2]=mod[1];
        mod[5]=mod[4];
        mod[8]=mod[7];
        break;
    case 406:   // Tetragonal c11 c33 c23 c12 c44 c66
        mod[1]=param[1];
        mod[2]=mod[1];
        mod[3]=param[2];
        mod[4]=param[3];
        mod[5]=mod[4];
        mod[6]=param[4];
        mod[7]=param[5];
        mod[8]=mod[7];
        mod[9]=param[6];
         break;
    case 407:   // Tetragonal c11 c33 c23 c12 c44 c66 c16
        mod[1]=param[1];
        mod[2]=mod[1];
        mod[3]=param[2];
        mod[4]=param[3];
        mod[5]=mod[4];
        mod[6]=param[4];
        mod[7]=param[5];
        mod[8]=mod[7];
        mod[9]=param[6];
        
        mod[12]=param[7];
         break;     
    case 306:   // Trigonal. 6 moduli
              // Give in order c11 c33 c23 c12 c44 c14        
         mod[1]=param[1];
         mod[2]=mod[1];
         mod[3]=param[2];
         mod[4]=param[3];
         mod[5]=mod[4];
         mod[6]=param[4];
         mod[7]=param[5];
         mod[8]=mod[7];
         mod[9]=0.5*(mod[1]-mod[6]);
         mod[10]=param[6];
        
         break;
    case 307:   // Trigonal. 7 moduli
              // Give in order c11 c33 c23 c12 c44 c14 c25
         mod[1]=param[1];
         mod[2]=mod[1];
         mod[3]=param[2];
         mod[4]=param[3];
         mod[5]=mod[4];
         mod[6]=param[4];
         mod[7]=param[5];
         mod[8]=mod[7];
         mod[9]=0.5*(mod[1]-mod[6]);
         mod[10]=param[6];
         mod[11]=param[7];
         
         break;
     
    case 9:   // Orthorhombic c11 c22 c33 c23 c13 c12 c44 c55 c66
         for(i=1;i<10;i++) mod[i]=param[i];
         
         break;   
         
    case 13:  // Monoclinic c11 c22 c33 c23 c13 c12 c44 c55 c66 c15 c25 c35 c46 
         for(i=1;i<14;i++) mod[i]=param[i];
    
    case 21:  // Triclinic c11 c22 c33 c23 c12 c44 c55 c66 c14 c15 c16 c24 c25 c26 c34 c35 c36 c45 c46 c56  
         for(i=1;i<22;i++) mod[i]=param[i];
             
         break;   
   
   }
   // Nye, J. F. Physical Properties of Crystals. London: Oxford University Press, 1959.
      c[1][1] =  10.*mod[1];
      c[1][2] =  10.*mod[6];
      c[1][3] =  10.*mod[5];
      
      c[1][4] =  10.*mod[10]; 
      c[1][5] = -10.*mod[11]; //Notice the - sign. The c25 is the input.
      c[1][6] =  10.*mod[12];
      
      c[2][1] = c[1][2];
      c[2][2] = 10.*mod[2];
      c[2][3] = 10.*mod[4];
      
      c[2][4] = -c[1][4]; 
      c[2][5] = -c[1][5];
      c[2][6] = -c[1][6];
      
      c[3][1] = c[1][3];
      c[3][2] = c[2][3];
      c[3][3] = 10.*mod[3];
      
      c[3][4] = 0;
      c[3][5] = 0;
      c[3][6] = 0;
      
      c[4][1] = c[1][4];
      c[4][2] = c[2][4]; 
      c[4][3] = 0;
      
      c[4][4] = 10.*mod[7];
      c[4][5] = 0;
      c[4][6] = c[2][5];
      
      c[5][1] = c[1][5];
      c[5][2] = c[2][5];
      c[5][3] = 0;
      
      c[5][4] = 0; 
      c[5][5] = 10.*mod[8];
      c[5][6] = c[1][4]; 
      
      c[6][1] = c[1][6];
      c[6][2] = c[2][6];
      c[6][3] = 0;
      
      c[6][4] = c[2][5];
      c[6][5] = c[1][4]; 
      c[6][6] = 10.*mod[9];
      
   if(iparm[36]==13) // Monoclinic diad along x2
   {
      c[1][4] =  0.; 
      c[1][5] = 10.*mod[10];
      c[1][6] =  0.;
      
      c[2][4] = 0.; 
      c[2][5] = 10.*mod[11];
      c[2][6] = 0.;
      
      c[3][4] = 0;
      c[3][5] = 10.*mod[12];
      c[3][6] = 0;
      
      c[4][1] = c[1][4];
      c[4][2] = c[2][4]; 
      c[4][3] = c[3][4];
      
      c[4][5] = 0;
      c[4][6] = 10.*mod[13];
      
      c[5][1] = c[1][5];
      c[5][2] = c[2][5];
      c[5][3] = c[3][5];
      
      c[5][4] = 0; 
      c[5][6] = 0; 
      
      c[6][1] = c[1][6];
      c[6][2] = c[2][6];
      c[6][3] = c[3][6];
      
      c[6][4] = c[4][6];
      c[6][5] = c[5][6]; 
   }
   
   
   if(iparm[36]==21) // Triclinic 
   {
      c[1][4] = 10.*mod[10]; 
      c[1][5] = 10.*mod[11];
      c[1][6] = 10.*mod[12];
      
      c[2][4] = 10.*mod[13]; 
      c[2][5] = 10.*mod[14];
      c[2][6] = 10.*mod[15];
      
      c[3][4] = 10.*mod[16];
      c[3][5] = 10.*mod[17];
      c[3][6] = 10.*mod[18];
      
      c[4][1] = c[1][4];
      c[4][2] = c[2][4]; 
      c[4][3] = c[3][4];
      
      c[4][5] = 10.*mod[19];
      c[4][6] = 10.*mod[20];
      
      c[5][1] = c[1][5];
      c[5][2] = c[2][5];
      c[5][3] = c[3][5];
      
      c[5][4] = c[4][5]; 
      c[5][6] = 10.*mod[21]; 
      
      c[6][1] = c[1][6];
      c[6][2] = c[2][6];
      c[6][3] = c[3][6];
      
      c[6][4] = c[4][6];
      c[6][5] = c[5][6]; 
   }
   
   
   
   
   
   free_vector(mod);
}

/* Get the frequencies from eigenvalues (divide by 2pi, etc. ) */ 
void get_freq (double *d, int dim)
{
int i;
for (i = 1; i <= dim; i++)
   {
   if (d[i] > 0) {d[i] = (sqrt (d[i])) / (2. * PI * 0.1 );}
   else d[i] = 0;
   }
}

void solve_freq(long *lookup,double *param, int *iparm, double *dfact,double **c,long *freq,long *fsym)
{
double **EK;
double **EP;
double *d;  // store diagonal elements 
double *e;  // store off-diagonal elements  
dimensions dimsample;
dimsample.dx=param[30];
dimsample.dy=param[31];
dimsample.dz=param[32];
dimsample.geom=iparm[45];

if(dimsample.geom==4)
    {
    dimsample.dxm=param[50];
    dimsample.dym=param[51];
    dimsample.dzm=param[52];    
    }

angles euler;
euler.a=remainder(param[33],360.0);
euler.b=remainder(param[34],360.0);
euler.c=remainder(param[35],360.0);
int i,j,dim,pos;
long RA[3];
double      *frth = vector(iparm[43]);   // Stores frequency
int         *ford = ivector(iparm[43]);   // Stores order w/in block
int         *frsym = ivector(iparm[43]);  // Stores block
index_frth *testf = indfrth(iparm[43]); //Used just for qsort.
// double time_spent=0.0;
// qsort(wordFreqs, numWords, sizeof(struct WordFrequency), wfcmp); 

pos=0;
i=1;
while(i<9)  // This could be parallelized - but only if symmety elements are present. Possibly parallel version of eigenvalue exists in LAPACK
   {
//    clock_t begin0 = clock();
    dim = ((int *)lookup[2])[i];
//    printf("%i %i\n",i,dim);
    gen_ke_pe_mat(RA,lookup,iparm[42],param[39],dimsample,c,i,euler,dfact); // 1
//    clock_t end0 = clock();
//    time_spent = (double)(end0-begin0)/CLOCKS_PER_SEC;
//    printf("Gen_KE_PE %f seconds |",time_spent);    
    
    EK = (double **) RA[0];
    EP = (double **) RA[1];
    
    d=vector(dim+1);
    e=vector(dim+1);
    
   // Plan to solve A.x = gB.x where A = EP and B = EK
   //   i)   calculate L and L' from B; L.L' = B
   //   ii)  calculate C = L^-1.A.(L^-1)'
   //   iii) calculate eigenvalues of C
   double **L = cholesky (EK, dim);
  // printf("i = %i  dim = %i\n",i,dim);
   
   
//   clock_t begin1 = clock();
   make_c (EP, L, dim); // 2 Now EP == C 
//   clock_t end1 = clock();
//   time_spent = (double)(end1-begin1)/CLOCKS_PER_SEC;
//   printf("Make c %f seconds |",time_spent);
   
//   clock_t begin2 = clock();
   householder (EP, dim, d, e); // 3 Tridiagonalize EP 
//   clock_t end2 = clock();
//   time_spent = (double)(end2-begin2)/CLOCKS_PER_SEC;
//   printf("Make householder %f seconds |",time_spent);
   
//   clock_t begin3 = clock();
   tri_eig (d, e, dim); // 4 Get the eigenvalues in e 
//   clock_t end3 = clock();
//   time_spent = (double)(end3-begin3)/CLOCKS_PER_SEC;
//   printf("Make tri_eig %f seconds | Dim = %i\n",time_spent,dim);
   
   
   for(j=1;j<=dim;j++)
      {
      frth[j+pos]=d[j];
      testf[j+pos].fr=d[j];
//      testf[j+pos].order=j;
      testf[j+pos].symm=i;
      //printf("%lg %lg %i %i\n",frth[j+pos],testf[j+pos].fr,testf[j+pos].order,testf[j+pos].symm);
      }
   pos=pos+dim;
   free_matrix (EP);
   free_matrix (EK);
   free_matrix (L);
   free_vector (d);
   free_vector (e);
   if(iparm[46]==0)  // If no mirror planes do only i=1. Else do 8 reps.If a rep has size 0, it will be skipped
       {i=9;}
       else
       {i++;}
   }
//ruserror("Stop before sort");
//qsort (frth,iparm[23], sizeof (double), compare_doubles);
qsort (testf,iparm[43], sizeof (index_frth), compare_indfrth);

for(i=1;i<iparm[43]+1;i++) {
    frth[i]=testf[i].fr;
// Add indexing here.
    ford[i]=testf[i].order;
    frsym[i]=testf[i].symm;
}

get_freq (frth,iparm[43]);
for(i=1;i<iparm[43]-5;i++) {
    frth[i]=frth[i+6];
    ford[i]=ford[i+6];
    frsym[i]=frsym[i+6];
} // Skip first 6 lines which are 0 for translations and rotations

freq[0] = (long) frth;
fsym[0] = (long) frsym;
free_indfrth(testf);
free_ivector(ford);

}

/*Get derivative of frequencies for LM routine */
void deriv_freq(long *lookup,double *param, int *iparm, double *dfact, double **c,long *freq,double *deltaparam,double **dfreq,int ideriv) {
int i,j;
double **cdelta = matrix (6);
double *dparam=vector(65);
long tempderiv[1];
long tempfsym[1];
// In all case ideriv should be <=15 and never be 12
for(i=1;i<65;i++) {dparam[i]=param[i];}
dparam[ideriv]+=deltaparam[ideriv];

if(ideriv==31||ideriv==30) {dparam[32]=param[32]*param[ideriv]/dparam[ideriv];} // Keep volume constant
                                                                                // Needs revision!
                                                       
gen_elatensor(dparam,iparm,cdelta);
solve_freq(lookup,dparam,iparm,dfact,cdelta,tempderiv,tempfsym);

                       for(j=1;j<=iparm[44];j++) 
                       {
                       ((double *)tempderiv[0])[j]-=((double *)freq[0])[j];
         		       ((double *)tempderiv[0])[j]/=deltaparam[ideriv];
          		       dfreq[ideriv][j]=((double *)tempderiv[0])[j];
		               }
free_matrix(cdelta);
free_vector(dparam);
free_ivector((int *) tempfsym[0]);
free_vector((double *) tempderiv[0]);
}

/* Get dependence of resonances on specific moduli */
void gen_dfdc (long *dfdc,long *lookup, double *param, int *iparm,double *dfact, long *fth,long *fsym)  
{
int i,j;
int nlines;
int nmods = iparm[36]%100;

if(iparm[38]==0) {nlines=iparm[44];}
    else {nlines= iparm[38];}

double buffer;
double *tempdfdc = vector (nlines*nmods+1);  
double epsilon = 0.001;
long fth_temp[1];
long fsym_temp[1];

printf("---- Generating derivatives.\n");
for (i = 1; i < nmods +1; i++)                        // This could be parallelized.
    {
    buffer=param[i];
    param[i]=buffer*(1.+epsilon);
    calclines (lookup,param,iparm,dfact,fth_temp,fsym_temp);            
    for(j = 1; j < nlines+1;j++)
        {
        tempdfdc[(i-1)*nlines+j] = 2.*( ((double *)fth_temp[0])[j]/((double *)fth[0])[j] - 1 )/epsilon;
        }
    param[i]=buffer;    
    
    dfdc[0]= (long) tempdfdc;
    free_vector((double *)fth_temp[0]);
    free_ivector((int *)fsym_temp[0]);
    }


}

double fchisq(long *fex, long *fth, long *w,long *fsym,int npts,int nfree) // w[0] is the list of weights.
    {
     int i;
     double chisq=0;
     if(npts==0||nfree==0) return 0.;     
     for(i=1;i<=npts;i++) 
        {
 	 chisq=chisq+((double *)w[0])[i]*SQR((((double *)fex[0])[i]/((double *)fth[0])[i])-1.);
//     printf("%g %g\n",((double *)fex[0])[i],((double *)fth[0])[i]);
	}
     return chisq/(0.0001*nfree); // Convert to percent.
	}

// Core fitting routine using Levenberg-Maqrquardt
int curfit (long *lookup, double *param, int *iparm, double *dfact,
             long *fex, long *fth, long *w, long *fsym,
	    double *chisqr0, double *chisqr1, double *flambda,  double *deltap,  double *sigmap  )
{
int i,j,k;
int ndata=0;
int nterms=0;
int nfree=0;
int ideriv,next;

int *dvar=ivector(35);
double **deriv = table(35,iparm[44]);
double *b=vector(66);
double **c = matrix (6); // Elastic Pseudo-Tensor

long tempfsym[1],tempfth[1];

j=1;
// Find the number of degrees of freedom of the problem
// Setup pairing of free parameters with derivative variable in dvar
// Counts the number of data points with w!=0

for(i=1;i<iparm[44]+1;i++)
   {
    if(  ((double *)w[0])[i]!=0.) ndata++; 
   }

nterms=iparm[36]%100; // Number of free parameters is number of elastic constants. Trig., tetra 306,307,406,407. 

for(i=1;i<iparm[36]%100+1;i++) 
{
    if( iparm[i]!=0) {nterms--;} // If an elastic constants is fixed, decrease by one
    else {dvar[j]=i;j++;}        // If it is free pair parameter i with dvar[j]
   }   
if(iparm[30]!=1){nterms+=2;dvar[j]=30;dvar[j+1]=31;j=j+2;} // If dimensions are free, 2 additional parameters. 12 stays constrained by 10 and 11.
if(iparm[33]!=1){nterms+=1;dvar[j]=33;j++;} // Euler angle psi
if(iparm[34]!=1){nterms+=1;dvar[j]=34;j++;} // Euler angle theta
if(iparm[35]!=1){nterms+=1;dvar[j]=35;j++;} // Euler angle phi
//print_ivector(dvar,nterms);
nfree=ndata-nterms; // real number of degrees of freedom
nfree=ndata; // similar to Migliori et al

iparm[41]=nterms;

//printf("Nfree %i Ndata %i Nterms %i\n",nfree,ndata,nterms);
if(nfree<=0||nterms==0) {ruserror("Nothing to fit");} 

// //////////////////////////////////////////////////////

// Get f(x), the theoretical frequencies.
gen_elatensor(param,iparm,c); 
solve_freq(lookup,param,iparm,dfact,c,fth,fsym); 

//Debugging: uncomment to see fth and fex //print_vector((double *)fth[0],iparm[24]);//print_vector((double *)fex[0],iparm[24]);
    
if(iparm[38]==0) { // If we are fitting

// Evaluate chisquare at starting point
*chisqr0=fchisq(fex,fth,w,fsym,iparm[44],nfree); // iparm[44] is the number of data points. 

// //////////////////////////////////////////////////////
printf("---- Initial Chisqr: %lg\n",*chisqr0);
//print_vector((double *)fex[0],iparm[44]);print_vector((double *)fth[0],iparm[44]);

double **array=matrix(nterms);
double **alpha=matrix(nterms);
double *beta=vector(nterms);

// Get alpha and beta 
    // dfreq[i] contains derivative  
	for(j=1;j<=nterms;j++)
           {
           ideriv=dvar[j]; // parameter w/r to which to get the derivative
           deriv_freq(lookup,param,iparm,dfact,c,fth,deltap,deriv,ideriv);   // This could be parallelized. Needs rewriting next loop also.
           for(i=1;i<=iparm[44];i++) 
                {
                beta[j]=beta[j]+ ((double *)w[0])[i]*( ((double *)fex[0])[i] - ((double *)fth[0])[i] ) * ((double *)deriv[ideriv])[i];	    
                for(k=1;k<=j;k++) alpha[j][k]=alpha[j][k]+(  ((double *)deriv[ideriv])[i] * ((double *)deriv[dvar[k]])[i] * ((double *)w[0])[i] );
	            }   
	       }
    free_vector((double *)fth[0]);
    free_ivector((int *)fsym[0]);
// Symetrize alpha    
for(j=2;j<=nterms;j++) for(k=1;k<j;k++) alpha[k][j]=alpha[j][k];
// //////////////////////////////////////////////////////    

// Invert modified curvature matrix to find new parameters
next=1;
while(next)
    {
    next=0;
    for(j=1;j<=nterms;j++)
		{
			for(k=1;k<=nterms;k++) {array[j][k]=alpha[j][k]/sqrt(alpha[j][j]*alpha[k][k]);}
			array[j][j]=array[j][j]*(1.+*flambda); //Modified.
		}


    // Invert array with Shipley-Coleman algorithm
    for(k=1;k<=nterms;k++)                
        {
        if(array[k][k]==0) {printf ("Error diagonal element %i=(%f)\n",k,array[k][k]);return 0;}
        array[k][k]= -1./(array[k][k]); 
        for(i=1;i<=nterms;i++) if (i!=k) { array[i][k]*=array[k][k]; } 
        for(i=1;i<=nterms;i++) if (i!=k) {for(j=1;j<=nterms;j++) if(j!=k) {array[i][j]+=(array[i][k])*(array[k][j]);}}
        for(i=1;i<=nterms;i++) if (i!=k) { array[k][i]*=array[k][k]; } 
        } 
	for(i=1;i<=nterms;i++) for(j=1;j<=nterms;j++) if (array[i][j]!=0.) {array[i][j]=-array[i][j];}
	// Inversion stops here. 
	
    for(i=1;i<65;i++) {b[i]=param[i];} // Use old parameters as starting value.
    for(j=1;j<=nterms;j++) for(k=1;k<=nterms;k++) b[dvar[j]]=b[dvar[j]]+beta[k]*array[j][k]/sqrt(alpha[j][j]*alpha[k][k]); // Get new solution
	
    // Special constrain for volume. Total volume/density needs to be conserved.
    b[32]=param[32]* (param[30]*param[31])/(b[30]*b[31]);
    if(iparm[45]==8)  b[32]=param[32]* (pow(param[30],2)-pow(param[31],2))/(pow(b[30],2)-pow(b[31],2)); // Hollow cylinder
    
    param[37]=lengthgrad(b,param,32);
    param[38]=*flambda;
    gen_elatensor(b,iparm,c);
    solve_freq(lookup,b,iparm,dfact,c,fth,fsym); // Get new frequencies for new solution.

    *chisqr1=fchisq(fex,fth,w,fsym,iparm[44],nfree);

    //printf("Initial chisq %f -> last obtained chisq %f\n",*chisqr0,*chisqr1);
    // If chi2 increased, increase flambda and try again
    if((*chisqr0-*chisqr1)<=-param[60]) {*flambda=*flambda*10.;next=1;printf("---- Step Chisqr:    %lg\n",*chisqr1);}
    }
	//printf("\n");
	
printf("---- Final Chisqr:   %lg\n",*chisqr1);


// Slow and correct error calculation
if(iparm[61]>0)
{	
	double buff=0,chi1=0,chi3=0,dp=0,chi2=*chisqr1;
	for(i=1;i<=nterms;i++)
		{
		buff=b[i];
		dp=10*(b[i]-param[i]);
		if(dp<1e-6){dp=1e-6;} // Make sure that parameter is changed by minimum amount.
		b[i]=buff+dp;
		calclines (lookup,b,iparm,dfact,tempfth,tempfsym);
		chi1=fchisq(fex, tempfth, w,tempfsym,iparm[44],nfree);
		free_vector((double *) tempfth[0]);
		free_ivector((int *) tempfsym[0]);
		b[i]=buff-dp;
		calclines (lookup,b,iparm,dfact,tempfth,tempfsym);
		chi3=fchisq(fex, tempfth, w,tempfsym,iparm[44],nfree);
		free_vector((double *) tempfth[0]);
		free_ivector((int *) tempfsym[0]);
		b[i]=buff;
		sigmap[i]=(fabs(dp)*2.*sqrt(0.01*chi2/(chi1-2*chi2+chi3)));
		}
}
else
{
	// Fast error estimation from curvature matrix.
	for(i=1;i<=nterms;i++) sigmap[i]= sqrt(array[i][i]/alpha[i][i])/10000.; 
}
	
for(i=1;i<65;i++) {param[i]=b[i]; if(b[i]<0) printf("==== Warning: parameter %i found negative.\n",i);} // New solution found. c

*flambda=*flambda/10.;
free_matrix(array);
free_matrix(alpha);
free_vector(beta);
      }
else
      {
      print_vector((double *)fth[0],iparm[38]); // Used only in the case that all parameters are fixed but a fit was attempted.
      *chisqr1=1.;
      }
        free_ivector(dvar);
        free_table(deriv);
        free_vector(b);
        free_matrix(c);
return 1;
}

// Core RUS frequency calculation.
int calclines (long *lookup, double *param, int *iparm,double *dfact,long *fth,long *fsym)
{
double **c = matrix (6); // Elastic Pseudo-Tensor
// Get f(x), the theoretical frequencies.
gen_elatensor(param,iparm,c); 
solve_freq(lookup,param,iparm,dfact,c,fth,fsym); 
free_matrix(c);
return 1;
}
