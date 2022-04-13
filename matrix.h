/* Declaration of functions used in matrix.c The puropse of each
   function is given in matrix.c file.*/

#ifndef _MATRIX_H
#define _MATH_H 1
#define PI  3.1415926535897932


// Definition of structures.
typedef struct
{
  double dx, dy, dz;
  double dxm, dym, dzm;
  double ID, OD;
  double logdx, logdy, logdz;
  int geom; // RPR by default
} dimensions;

typedef struct
{
  double   fr;
  int      order,symm;
} index_frth;

typedef struct
{
  int   block[32];
  char  comment[32][128];
} comments;

typedef struct
{
  double chi2;
  double param_c[32];
  double delta_c[32]; 
} solution;

typedef struct 
{
  int i, l, m, n;
} indices;

typedef struct 
{
  double a,b,c;
} angles;


// Mathematical functions.
double delta  (int, int);
double integ (int xexp, int yexp, int zexp, dimensions a,double *dfact);
double pythag (double a, double b);

// IO stream and error functions
void ruserror (char *);
void println (void);
void print_matrix (double **mat, int dim);
void print_vector (double *vec, int dim);
void print_vector_h (double *vec, int dim);
void print_ivector (int *vec, int dim);
void print_output (double *vec, int no,int n);

void read_input (char *filename, comments *COMM, int *iparm, double *param,long *fex,long *w);
void write_io (char *filename, comments *COMM, int *iparm, double *param,long *fex,long *fth,long *w,long *fsym);
void write_out (char *filename, comments *COMM, int *iparm, double *dfact, double *param,long *fex,long *fth, long *w,long *fsym,double *sigmap, long *dfdc, long *lookup);

// Comparison functions for qsort.
int compare_doubles (const void *a, const void *b); 
int compare_indfrth (const void *f1, const void *f2); 
int compare_solution (const void *f1, const void *f2); 

// Object definition, memory allocation and memory cleaning functions.
double **matrix (int dim);
void free_matrix (double **mat);

index_frth *indfrth(int dim);
void free_indfrth(index_frth *indfrth);

solution *sol(int dim);
void free_solution(solution *sol);

double     *vector (int dim);
void free_vector (double *vec);

int *ivector (int dim);
void free_ivector (int *vec);

// Matrix functions.
void mul_matrix_sca (double **mat, double sca, int dim);
double *mul_matrix_vec (double **mat, double *vec, int dim);
double **mul_matrix_matrix (double **mata, double **matb, int dim);
double **transpose_matrix (double **mat, int dim);
void householder (double **a, int n, double *d, double *e);
void tri_eig (double *d, double *e, int n);
double **cholesky (double **mat, int dim);
void make_c (double **a, double **l, int dim);

// RUS specific functions
void gen_dfdc (long *dfdc,long *lookup, double *param, int *iparm,double *dfact,long *fth,long *fsym);
void get_freq (double *d, int dim);
void gen_ke_pe_mat(long *RA,long *lookup,int nmax,double rho,dimensions dimsample,double **c,int bdi, angles euler, double *dfact);
void gen_lookups(long *lookup,int nmax, int mirr);
void gen_elatensor(double *param, int *iparm,double **c);
void solve_freq(long *lookup,double *param, int *iparm,double *dfact, double **c,long *freq,long *fsym);
void deriv_freq(long *lookup,double *param, int *iparm,double *dfact, double **c,long *freq,double *deltaparam,double **dfreq, int ideriv);
int calclines (long *lookup, double *param, int *iparm,double *dfact, long *fth,long *fsym);	   
indices getind (int, int);


// Levenberg-Marquardt functions.
int curfit (long *lookup, double *param, int *iparm,double *dfact, long *fex, long *fth, long *w,long *fsym,double *chisqr0, double *chisqr1 ,
	   double *flambda, double *deltap, double *sigmap);
double fchisq(long *fex, long *fth, long *w,long *fsym,int npts,int nfree);
double lengthgrad(double *vec, double *vec2, int dim);
	   
#endif /* math.h */ 
