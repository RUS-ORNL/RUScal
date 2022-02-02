/* Declaration of functions used in matrix.c The puropse of each
   function is given in matrix.c file.*/

#ifndef _MATRIX_H
#define _MATH_H 1

#define PI  3.1415926535897932

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
  int i, l, m, n;
} indices;


typedef struct 
{
  double a,b,c;
} angles;

void ruserror (char *);
void println (void);

indices getind (int, int);
double delta  (int, int);
double integ (int xexp, int yexp, int zexp, dimensions a,double *dfact);

double fchisq(long *fex, long *fth, long *w,long *fsym,int npts,int nfree);

double pythag (double a, double b);

int compare_doubles (const void *a, const void *b); 
int compare_indfrth (const void *f1, const void *f2); 


double **matrix (int dim);
void free_matrix (double **mat);

index_frth *indfrth(int dim);
void free_indfrth(index_frth *indfrth);

double     *vector (int dim);
void free_vector (double *vec);

int *ivector (int dim);
void free_ivector (int *vec);

void print_matrix (double **mat, int dim);
void print_vector (double *vec, int dim);
void print_vector_h (double *vec, int dim);
void print_ivector (int *vec, int dim);
void print_output (double *vec, int no,int n);
double lengthgrad(double *vec, double *vec2, int dim);

void mul_matrix_sca (double **mat, double sca, int dim);
double *mul_matrix_vec (double **mat, double *vec, int dim);
double **mul_matrix_matrix (double **mata, double **matb, int dim);
double **transpose_matrix (double **mat, int dim);

void householder (double **a, int n, double *d, double *e);
void tri_eig (double *d, double *e, int n);
double **cholesky (double **mat, int dim);
void make_c (double **a, double **l, int dim);

//void read_input (char *filename, int *nmax, double *rho, double *dx,double *dy, double *dz, double *modulus,long *fex,long *wt,int *n,int *nummod);
void read_input (char *filename, comments *COMM, int *iparm, double *param,long *fex,long *w);
void write_io (char *filename, comments *COMM, int *iparm, double *param,long *fex,long *fth,long *w,long *fsym);
void write_out (char *filename, comments *COMM, int *iparm, double *dfact, double *param,long *fex,long *fth, long *w,long *fsym,double *sigmap, long *dfdc, long *lookup);
void gen_dfdc (long *dfdc,long *lookup, double *param, int *iparm,double *dfact,long *fth,long *fsym);

void gen_ke_pe_mat(long *RA,long *lookup,int nmax,double rho,dimensions dimsample,double **c,int bdi, angles euler, double *dfact);
void gen_lookups(long *lookup,int nmax, int mirr);
void gen_elatensor(double *param, int *iparm,double **c);
void solve_freq(long *lookup,double *param, int *iparm,double *dfact, double **c,long *freq,long *fsym);
void deriv_freq(long *lookup,double *param, int *iparm,double *dfact, double **c,long *freq,double *deltaparam,double **dfreq, int ideriv);

void get_freq (double *d, int dim);

int curfit (long *lookup, double *param, int *iparm,double *dfact, long *fex, long *fth, long *w,long *fsym,double *chisqr0, double *chisqr1 ,
	   double *flambda, double *deltap, double *sigmap);

int calclines (long *lookup, double *param, int *iparm,double *dfact, long *fth,long *fsym);	   
	   
#endif /* math.h */ 
