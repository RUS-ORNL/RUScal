/* rus.c: main file for RUScal -- a program to calculate resonant frequencies
*  of crystals and polycrystals
* with rpr, cylinder, ellipsoid, octahedral, or hollow cylinder geometry
* with rotations
* with or without mirror planes
* tools include simulation, fit, Monte-Carlo, genetic algorithm and grid search.
* main code does not contain any specifc geometry or symmetry information
* all calculations are done in curfit() or in calclines() - see matrix.c
* all functions are in matrix.c - main is only once-through.
*/ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "matrix.h"

int main (int argc, char *argv[])
{
    // First argument is the input filename if none, use 'rusin.dat'
    char *filename;
	if (argc < 2) {filename="rusin.dat";} else {filename = argv[1];}

	// Prepare memory for parameters and constraints
    double *param = vector(65);
    int    *iparm = ivector(65);
    double *deltap = vector(65);
    double *sigmap = vector(65);
	
    double *dfact  = vector(146); // Store double factorial values. Shifted by 2. dfact[1] == (-1)!!; dfact[2] == 0!!; dfact[n+2] == n!! and factorial values in space 100-120
   
	double *mclc = vector(22); //Lower bound for MonteCarlo
    double *mcdc = vector(22); //Difference from lower to upper for MC
	
    comments COMM; // Stores all comments in input/output files.
    
    // Lookup table that maps polynomial (l,m,n) to sub-matrix
    long lookup[3];
    long fex[1],fth[1],fsym[1],w[1],dfdc[1]; // Pointer to list of exp, th, symm, weight and derivative for frequencies.

    // Controls fo chi-square minization routine
    double flambda=0.0001;           // Control of iterations  
    double chisqr1=1.0;           // Chisquare - will contain the chisqr after calls to curfit
    double chisqr0=1.0;           // Chisquare - will contain the initial chisqr after calls to curfit
    
	int i,j,k,l,cd;  // Helper variables
	double cbuff,sum,b_mass,b_param; // Helper variables
    
    printf("\n======================================\n");	
    printf("Please cite J. Torres et al., JASA ?? (202?)\n");
    printf("======================================\n");	  
	
	read_input (filename,&COMM,iparm, param,fex,w);
	
    if(iparm[42]>20) {ruserror("Polynomial can't be larger than 20.\n");}
    else
    {   // Generate double factorials and factorials. dfact[n+2]=n!! ; dfact[n+102]=n!
        dfact[0]=1; // Sanity    
        dfact[1]=1; // (-1)!!
        dfact[2]=1; // (0)!!
        dfact[3]=1; // (1)!!
        for(i=4;i<2*iparm[42]+6;i++) {dfact[i]=dfact[i-2]*(i-2);}
        dfact[98]=1;
        dfact[99]=1;
        dfact[100]=1;
        dfact[101]=1;
        for(i=2;i<2*iparm[42]+6;i++) {dfact[100+i]=dfact[100+i-1]*i;}
    }
    
    gen_lookups(lookup,iparm[42],iparm[46]); // Generate lookup table

	// By default initial variation of parameter variation is 1/1000.
	for(i=1;i<36;i++) 
	{
		deltap[i]=param[i]*0.001; // Standard for initial steps
		if(param[i]==0) {deltap[i]=0.01;} // Should happen only if fitting Euler angle with initial 0 value. Choose abritrary positive step.
	} 

	printf("\n======================================\n");	      

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initializations are finished.
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Following action depends on chosen mode.

   	if(iparm[38]>0)   // Simple calculation. Calculate lines, number given by iparm[38].
	{
		FILE *fp;
		if ((fp = fopen ("freqout_calc.dat", "a")) == NULL) ruserror ("cannot open output file");
		fprintf(fp,"%s",COMM.comment[0]);  // Write title line    
		
		calclines (lookup,param,iparm,dfact,fth,fsym);
		printf("==============================\n");
		
		if(iparm[39]) // If calculation of derivative of frequencies vs elastic constants was requested.
		{
			gen_dfdc(dfdc,lookup,param,iparm,dfact,fth,fsym); 
			for(j=1;j<iparm[38]+1;j++)
			{
				printf("%3i %8.6f",j,((double *)fth[0])[j]); 
				fprintf(fp,"%3i %8.6f",j,((double *)fth[0])[j]); 
				sum=0;
				for(i=1;i<iparm[36]%10+1;i++)
				{
					sum+= ((double *)dfdc[0])[(i-1)*iparm[38]+j];  
					printf("%5.2f ",((double *)dfdc[0])[(i-1)*iparm[38]+j]);
					fprintf(fp,"%5.2f ",((double *)dfdc[0])[(i-1)*iparm[38]+j]);
				}
				printf("%5.2f\n",sum);
				fprintf(fp,"%5.2f\n",sum);
			}
			free_vector((double *)dfdc[0]);
		}
		else
		{
			for(j=1;j<iparm[38]+1;j++)
			{
				printf("%3i %8.6f\n",j,((double *)fth[0])[j]); 
				fprintf(fp,"%3i %8.6f\n",j,((double *)fth[0])[j]); 
			}
		}
	fclose(fp);
	free_vector((double *)fth[0]);
	free_ivector((int *)fsym[0]);        
	}
	
	if(iparm[38]==-1||iparm[38]==-5) // Monte carlo simulation (Mode -1) with optional genetic algorithm (Mode -5).
	{ 
		// Control parameters.
		int    sizeofmc=iparm[60]; // Number of MC steps
		int    numgen=iparm[61];    // Number of generations
		double selection=param[60]; // Environmental pressure
		double mut_range=param[61]; // Range of mutation. If 0: interpolate between paretns only. >0 allows extrapolation.
		double mut_prob=param[62]; // Mutation probability
		
		int nc=iparm[36]%100; // Number of elastic constants

		// Helper variables
		FILE *fp;
		double chi2;
		int ndata=0;
		int elite=0;
		int i1,i2;
		double breed,mut;
		
		// Initialize lists that will serve to store populations.
		solution *sol1 = sol(sizeofmc);
		solution *pop = sol(sizeofmc);
		
		printf("#= Starting Monte-Carlo==========================\n");
		
		for(i=1;i<iparm[44]+1;i++)
		{
			if(  ((double *)w[0])[i]!=0.) ndata++; // Count degrees of freedom
		}

		time_t t;
		srand((unsigned) time(&t)); // Initialize random number generation with time as seed.

		if ((fp = fopen ("rms_out.dat", "a")) == NULL) ruserror ("cannot open output file"); // Stores MC data on the fly.
		fprintf(fp,"%s",COMM.comment[0]);  // Write title line
		fprintf(fp,"#Step  rms     cxx\n");
		fprintf(fp,"#============================================\n");
		fclose(fp);

		for(i=1;i<=nc;i++)
		{
			mclc[i]=param[i]*(1.-iparm[i]/1000.); // Lower bound for parameter i
			mcdc[i]=param[i]*2.*iparm[i]/1000.;   // Delta for parameter i
		}
				
		for(j=0;j<sizeofmc;j=j+1)
		{
			for(i=1;i<=iparm[36]%10;i++)
			{
			param[i]=mclc[i]+mcdc[i]*(double)rand()/(double)RAND_MAX;
			sol1[j].param_c[i]=param[i]; // Store value of tested solution
			}
		
		calclines (lookup,param,iparm,dfact,fth,fsym);
		chi2=fchisq(fex, fth, w,fth,iparm[44],ndata);
		free_vector((double *)fth[0]); // Clean memory after calculating lines.
		free_ivector((int *)fsym[0]);  // Clean memory after calculating lines.
		
		sol1[j].chi2 = chi2; // Store chi2 for tested solution

		if ((fp = fopen ("rms_out.dat", "a")) == NULL) ruserror ("cannot open output file");
		fprintf(fp,"%4i %10.6f ",j,sqrt(chi2));
		for(i=1;i<=nc;i++)
			{
			fprintf(fp,"%8.6f ",param[i]); // Prints value of tested parameters.
			}
		fprintf(fp,"\n");
		fclose(fp);
		
		if(j%50==0)
			{
			printf(" %i\n",j);	// Progress bar on screen.
			}
		}

		qsort (sol1,sizeofmc, sizeof (solution), compare_solution); // Sort population by value of chi2.
		if ((fp = fopen ("rms_sort.dat", "a")) == NULL) ruserror ("cannot open output file");
		fprintf(fp,"%s",COMM.comment[0]);  // Write title line
		fprintf(fp,"#Step  rms     cxx\n");
		fprintf(fp,"#============================================\n");
		fprintf(fp,"#==========================\n");
		fprintf(fp,"#Sorted solutions at gen 0:\n");
		fprintf(fp,"#==========================\n");
		for(j=0;j<sizeofmc;j=j+1)
		{
			fprintf(fp,"%4i %10.6f ",j,sqrt(sol1[j].chi2));
			pop[j].chi2=sol1[j].chi2;
			for(i=1;i<=nc;i++)
			{
			fprintf(fp,"%8.6f ",sol1[j].param_c[i]);
			pop[j].param_c[i]=sol1[j].param_c[i]; // Use this loop to make a copy for the next generation.
			}
		fprintf(fp,"\n");
		}
		fclose(fp);
		// Now pop and sol1 contain the same info.
		
		if(iparm[38]==-5) 
		{
			printf("#= Starting Genetic algorithm ========================\n");
			printf("# Controls - Pop: %i Numgen: %i EvolP: %5.3f Extrapol: %5.3f Mut: %5.3f \n",sizeofmc,numgen,selection,mut_range,mut_prob);
		}
		else 
		{
			numgen=0; // Skip genetic algorithm.
		}
		
		for(l=0;l<numgen;l++) // Loop over generations
		{
			elite=0;	// A priori no elite parents are saved. Exception is made when all children perform worse.
			printf("=======================\n");
			printf("Starting generation %i:\n",l+1);
			for(j=0;j<sizeofmc;j++) // Generate new population
			{
				//Select individual based on three random numbers. 1->parent1, 2-> parent2, 3->breeding: how much of parent1 genes.
				i1=sizeofmc*(pow((double)rand()/(double)RAND_MAX,selection)); 
				i2=sizeofmc*(pow((double)rand()/(double)RAND_MAX,selection));
				breed=(double)rand()/(double)RAND_MAX;
				breed=breed*(1.+2.*mut_range)-mut_range; // Option: extend range from [0:1]  to [-0.x:1.x] (0<x<0.5)
			
				for(i=1;i<=nc;i++) // Loop over number of elastic constants.
				{
					mut=(double)rand()/(double)RAND_MAX; // Probability of mutation.
					if(mut<(mut_prob))		// Mutation frequency; default 1/number of constants. 
					{
						mut=( (mut*mut_prob)-0.5 )/(l+1); // Mutation range is -0.5<mut<0.5 of initial box (below iparm[i]/1000.) at step 0 then gets reduced at next generation.
					}  
				else
					{
						mut=0; // No mutation for this elastic constant.
					} 
				// Generate new solution.
				sol1[j].param_c[i]=pop[i1].param_c[i]*breed+pop[i2].param_c[i]*(1-breed)+mut*iparm[i]/1000.;
				param[i]=sol1[j].param_c[i]; // Prepares parameter vector for next calculation of lines.
				}
				calclines (lookup,param,iparm,dfact,fth,fsym);
				chi2=fchisq(fex, fth, w,fth,iparm[44],ndata);
				free_vector((double *)fth[0]); // Clean memory after calculating lines.
				free_ivector((int *)fsym[0]);  // Clean memory after calculating lines.

				sol1[j].chi2 = chi2; // Store chi2 for tested solution			
				if(j%50==0) 
				{
					printf(" %i\n",j); // Progress bar.
				}
			}
			qsort (sol1,sizeofmc, sizeof (solution), compare_solution);

			//Find cross-over point for elite parents.
			// Keep parents that are fitter than all children. Eliminate least fit individuals.
			while(sol1[0].chi2>pop[elite].chi2) 
			{
				elite++;
			}
			if(elite!=0)
			{
				k=0;
				for(j=sizeofmc-elite;j<sizeofmc;j++)
				{
					sol1[j].chi2=pop[k].chi2;
					for(i=1;i<=nc;i++)
					{
						sol1[j].param_c[i]=pop[k].param_c[i];
					}
					k++;
				}
				qsort (sol1,sizeofmc, sizeof (solution), compare_solution);	
			}

			
			if ((fp = fopen ("rms_sort.dat", "a")) == NULL) ruserror ("cannot open output file");
			fprintf(fp,"\n\n");
			fprintf(fp,"Best solutions at gen %i:\n",l+1);
			fprintf(fp,"==========================\n");
			for(j=0;j<sizeofmc;j=j+1)
			{
				fprintf(fp,"%4i %10.6f ",j,sqrt(sol1[j].chi2));
				pop[j].chi2=sol1[j].chi2;
				
				for(i=1;i<=nc;i++)
				{
					fprintf(fp,"%8.6f ",sol1[j].param_c[i]);
					pop[j].param_c[i]=sol1[j].param_c[i];
				}
				fprintf(fp,"\n");
			}
			fclose(fp);
		}

		// Genetic algorithm is done. Prepare for one final classical fit.
		for(i=1;i<=nc;i++)
		{
			param[i]=pop[0].param_c[i]; // Choose best solution.
			k=iparm[i];
			if(k!=0)
			{
				iparm[i]=0; // Parameters with mc interval become free
			}  
			else
			{
				iparm[i]=1; // Parameters with no mc interval become fixed
			}      
		}		
		iparm[38]=0; // Toggle back to fit mode.
		param[60]=0.00001; // Toggle chi2 stopping precision to 1E-5.
		iparm[61]=1; // Toggle error analysis to rigorous.
		free_solution(sol1); // Clean memory
		free_solution(pop);  // Clean memory
	}

	if(iparm[38]==-2) // Systematic variation in a dimension
	{ 
		FILE *fp,*seqmat; 
		double sum; 
		// Calculate lines number given by iparm[44], for a systematic variation of dimensions, at constant density.
		// The dimension to change is given by iparm[37]. 
		b_mass=param[40];
		// Control dimension is 30-32 (a length) or 50-52 (an Euler angle); set in cd.
		cd=0;
		iparm[38]=iparm[44]; // Number of lines in input file.
		if(iparm[37]>0&&iparm[37]<4)  cd=iparm[37]+29; // 1-3: x+,y+, or z+
		if(iparm[37]>3&&iparm[37]<7)  cd=iparm[37]+46; // 4-6: x-,y-, or z-
		if(iparm[37]>6&&iparm[37]<10) cd=iparm[37]+26; // 7-9: Euler angle.
		if(cd==0) ruserror("Specify which dimension (1-6) or angle (phi=7; th =8; psi = 9) to vary.\n");
		
		b_param=param[cd]; // Buffer initial parameter.
		if(b_param==0) ruserror("Variation of dimension or angle can't be done for a parameter equal to zero.\n");
		
		if ((fp = fopen ("sequence_out.dat", "w")) == NULL) ruserror ("cannot open output file");
		if ((seqmat = fopen ("sequence_matrix.dat", "w")) == NULL) ruserror ("cannot open output file");
		fprintf(seqmat,"%s",COMM.comment[0]);  // Write title line
		fprintf(fp,    "%s",COMM.comment[0]);  // Write title line
		fprintf(seqmat,"#Parameter|");
		for(i=1;i<=iparm[44];i++)  // Loop on number of lines to calculate. Header for seqmat file.
		{
			fprintf(seqmat,"Line %3i|",i);            
		}
		fprintf(seqmat,"\n");
	
		for(i=1;i<21;i++)
		{
			param[cd]=b_param/10.*i; // 20 steps from 10% to 200% of the initial value. For angles: advise is to use 100 deg, which calculates 10-200 deg.
			// Line below changes the mass if a dimension is changed. For Euler angles, no change.
			if(iparm[37]<7)
			{
				param[40]=b_mass/10.*i;  // For most shapes, mass is linear in the dimensions.
				switch(iparm[45])
				{
					// Use mass = density * volume   
					case 4: param[40] = param[39]/1000.*(param[30]+param[50])*(param[31]+param[51])*(param[32]+param[52])/(6/PI); break;
					case 8: param[40] = param[39]/1000.*(param[32]*PI/4*(param[30]*param[30]-param[31]*param[31])); break; // Hollow cylinder. OD**2 - ID**2.
				}
			}
						
			if(param[40]>0) // Unless mass is negative calculate lines.
			{
				calclines (lookup,param,iparm,dfact,fth,fsym);
				printf("=====================================\n");
				printf("= Dimension: %f, Mass: %f\n",param[cd],param[40]);
				printf("=====================================\n\n");
			
				print_vector((double *)fth[0],iparm[44]);
				if(iparm[39]) 
				{
					gen_dfdc(dfdc,lookup,param,iparm,dfact,fth,fsym); // Get derivative of frequencies vs moduli.
				}
				
				// Printout results in two files
				// One file per parameter for frequencies with dfdc; one block per parameter.
				// One summary file for frequencies vs control parameter.
				fprintf(fp,"# Param %i =%8.5f \n",cd,param[cd]);
				fprintf(seqmat," %8.5f ",param[cd]);
				for(j=1;j<=iparm[44];j++)  // Loop on number of lines to calculate.
				{
					sum=0; //added
					fprintf(fp, "%3i    %8.7f    ", j, ((double *)fth[0])[j] );
					fprintf(seqmat,"%8.5f ",((double *)fth[0])[j]);
					if(iparm[39])
						{
						for(k=1;k<iparm[36]%100+1;k++) // Loop on derivatives vs moduli.
							{
							sum+= ((double *)dfdc[0])[(k-1)*iparm[44]+j];
							fprintf(fp,"%8.5f    ",((double *)dfdc[0])[(k-1)*iparm[44]+j]);
							}
						fprintf(fp,"%8.5f",sum);
						}
					fprintf(fp,"\n");
				}
				fprintf(fp,"###                          ###\n");
				fprintf(fp,"\n");
				fprintf(seqmat,"\n");
				if(iparm[39]) {free_vector((double *)dfdc[0]);} // Clean memory
				free_vector((double *)fth[0]);
				free_ivector((int *)fsym[0]);
			}
			else 
			{
				// Do nothing except printing a warning. This happens only for a hollow cylinder if ID>OD.    
				fprintf(fp,"# ============================================= #\n");
				fprintf(fp,"# Skipping unphysical parameter value %8.5f . #\n",param[cd]);
				fprintf(fp,"# ============================================= #\n");
				fprintf(seqmat,"# Value %8.5f is unphysical\n",param[cd]);
			}
		}
		fclose(fp); 
		fclose(seqmat);
		iparm[38]=-2; // Reset to appropriate mode.
	}
	
	if(iparm[38]==-3) // Check out error propagation
	{ // Calculate lines
		FILE *fp;
		double dj=0;
		int nfree=0;
		if ((fp = fopen ("errorout.dat", "a")) == NULL) ruserror ("cannot open output file");
		for(i=0;i<iparm[44];i++) {	if( ((double *)w[0])[i]!=0.0) nfree++;	}
		fprintf(fp,"%s",COMM.comment[0]);  // Write title line
		fprintf(fp,"#  c_xx      chi2 (dc/c)_xx\n");
		for(j=1;j<=iparm[36]%100;j++)
		{	
			fprintf(fp,"# Deviation of elastic constant %i\n",j);
			printf(  "\n# Deviation of elastic constant %i\n",j);
			cbuff=param[j];
			for(i=-5;i<=5;i++)
			{
				if(i==0) {dj=0;} else {dj=i*i*i/abs(i)/10000.;} // Relative deviation is sign(i)*i^2/10000.
				param[j]=cbuff*(1.+dj);
				calclines (lookup,param,iparm,dfact,fth,fsym);				
				chisqr0=fchisq(fex, fth, w,fth,iparm[44],nfree);
				free_vector((double *)fth[0]);
				free_ivector((int *)fsym[0]);
				fprintf(fp,"%8.6f %8.6f % 8.6f\n",param[j],chisqr0,dj);
				printf("%8.6f %8.6f %8.6f\n",param[j],chisqr0,dj);
			}
			param[j]=cbuff;
			fprintf(fp,"\n\n");
		}
		printf("\n");
		fclose(fp);
		
	}

	if(iparm[38]==-4) // Grid calculation - up to three elastic constants
	{ // Calculate lines
		FILE *fp;
		double chi2;
		int ndata=0,ngrid=0,gs=1,imax,jmax,kmax,l;
		
		for(i=1;i<iparm[44]+1;i++)
		{
			if(  ((double *)w[0])[i]!=0.) ndata++; 
		}

		iparm[61]=0;iparm[62]=0;iparm[63]=0;
		for(i=1;i<=iparm[36]%100;i++)
		{
			if(iparm[i]!=0 && ngrid<3)    // Only the first three intervals for constants are used.
			{
				ngrid++;
				iparm[60+ngrid]=i;       // Store which parameters to grid in iparm[61,62,63]
			}
			else
			{
				iparm[i]=0;		
			}
			mclc[i]=param[i]*(1.-iparm[i]/1000.); // Lower bound for grid
			mcdc[i]=param[i]*2.*iparm[i]/1000.;   // Interval for grid
		}
		if (ngrid==0) ruserror ("no parameter to construct grid on");
		if ((fp = fopen ("chi2out_grid.dat", "a")) == NULL) ruserror ("cannot open output file");

		// Handle loops: which cxx to use for grid
		// Format output according to ngrid
		// ngrid directs how many nested loops to do.
		fprintf(fp,"%s",COMM.comment[0]);  // Write title line
		fprintf(fp,"#Step  chi2      c[%i]",iparm[61]);
		if (ngrid>1) fprintf(fp,"      c[%i]",iparm[62]);
		if (ngrid>2) fprintf(fp,"      c[%i]",iparm[63]);
		fprintf(fp,"\n");
		fprintf(fp,"#==============================================================\n");
		fclose(fp);

		gs=0; //gs = grid steps.
		imax=20*(iparm[61]!=0)+1;
		jmax=20*(iparm[62]!=0)+1;
		kmax=20*(iparm[63]!=0)+1;
		i=0;
		while(i<imax)
		{
			param[iparm[61]]=mclc[iparm[61]]+mcdc[iparm[61]]*i/20.;
			j=0;
			while(j<jmax)
			{
				param[iparm[62]]=mclc[iparm[62]]+mcdc[iparm[62]]*j/20.;
				k=0;
				while(k<kmax)
				{
					param[iparm[63]]=mclc[iparm[63]]+mcdc[iparm[63]]*k/20.;
					calclines (lookup,param,iparm,dfact,fth,fsym);
					chi2=fchisq(fex, fth, w,fth,iparm[44],ndata);
					free_vector((double *)fth[0]);
					free_ivector((int *)fsym[0]);
					if ((fp = fopen ("chi2out_grid.dat", "a")) == NULL) ruserror ("cannot open output file");
					fprintf(fp,"%5i %10.6f ",gs,chi2);
					for(l=1;l<=ngrid;l++)
					{
						fprintf(fp,"%8.6f ",param[iparm[60+l]]);
					}
					fprintf(fp,"\n");
					k++;
					gs++;
					fclose(fp);
				}   
				j++;
			}
			i++;
		}

	}
	
	if(iparm[38]==0)  // Fit
	{ 
		for(i=0;i<20;i++)
		{ // Limit to 20 iterations
			if(i<20)
			{
				printf("Iteration %i\n--------------------------------------\n",i);
				if (curfit (lookup,param,iparm,dfact,fex,fth,w,fsym,&chisqr0,&chisqr1,&flambda,deltap,sigmap) ) 
				{
					printf("---- %lg ",param[1]);
					for (j=2;j<iparm[36]%100+1;j++) printf("%lg ",param[j]);
					printf("\n");
				}
			}
			// Add logfile here with parameters and chisq	 
			//printf("=====================\nrms error = %lg%%\nChisq0 %lg Chisqr1 %lg Lambda %lg\n",sqrt(chisqr0),chisqr0,chisqr1,flambda);
			param[41]=sqrt(chisqr1);
			param[36]=1.-(chisqr1/chisqr0); // Store relative change in chisqr
			iparm[40]=i;
			printf("---- rms error = %lg%%\n",param[41]);
			printf("======================================\n");
			if(param[36]<param[60]) {i=19;} // When chisqr changes less than 1/10000 or set tolerance stop at next iteration. 
			write_io(filename,&COMM,iparm, param,fex,fth,w,fsym);        // Generate output files for next run
			chisqr0=chisqr1;
			write_out(filename,&COMM,iparm,dfact,param,fex,fth,w,fsym,sigmap,dfdc,lookup); // Generate output
			free_vector((double *)fth[0]);
			free_ivector((int *)fsym[0]);
		}
	}
	
	// Clean up memory.
	free_vector(param);
	free_vector(deltap);
	free_vector(sigmap);
	free_ivector(iparm);
	free_vector(mclc);
	free_vector(mcdc);
	free_vector(dfact);
	free_ivector((int *)lookup[0]);
	free_ivector((int *)lookup[1]);
	free_ivector((int *)lookup[2]);
	free_vector((double *)fex[0]);
	free_vector((double *)w[0]);
	printf("---- Run ended\n");
	printf("======================================\n\n");	
	return 1;
}
