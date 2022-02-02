/* rus.c: main file for RUScal -- a program to calculate resonant frequencies
*  of crystals and polycrystals
* with rpr, cylinder, ellipsoid, octahedral, or hollow cylinder geometry
* with rotations
* with or without mirror planes
* main code does not contain any specifc geometry or symmetry information
* all calculations are done in curfit() or in calclines() - see matrix.c
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
    comments COMM;

    if (argc < 2) {filename="rusin.dat";}
    else {filename = argv[1];}

    int i,j,k,cd;
    // Prepare memory for parameters and constraints
    double *param = vector(55);
    double *mclc = vector(22); //Lower bound for MonteCarlo
    double *mcdc = vector(22); //Difference from lower to upper for MC
    int    *iparm = ivector(65);
    double *deltap = vector(55);
    double *sigmap = vector(55);
    double *dfact  = vector(146); // Store double factorial values. Shifted by 2. dfact[1] == (-1)!!; dfact[2] == 0!!; dfact[n+2] == n!! and factorial values in space 100-120
   
    
    // Lookup table that maps polynomial (l,m,n) to sub-matrix
    long lookup[3];
    long dfdc[1];
    long fex[1],w[1],fth[1],fsym[1];

    // Controls fo chi-square minization routine
    double flambda;           // Control of iterations  
    double chisqr1;           // Chisquare - will contain the chisqr after calls to curfit
    double chisqr0;           // Chisquare - will contain the initial chisqr after calls to curfit
    flambda=0.0001;
    chisqr1=1.0;
    chisqr0=1.0;
    double cbuff;
    double sum;
    double b_mass,b_param;
    
    printf("\n======================================\n");	
    printf("Please cite J. Torres et al., JASA ?? (202?)\n");
    read_input (filename,&COMM,iparm, param,fex,w);
    printf("======================================\n");	      
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
    
    // Generate lookup table
    gen_lookups(lookup,iparm[42],iparm[46]); // Do always but only once

    // By default initial variation of parameter variation is 1/1000.
    for(i=1;i<36;i++) 
    {
        deltap[i]=param[i]*0.001; // Standart for initial steps
        if(param[i]==0) {deltap[i]=0.01;} // Should happen only if fitting Euler angle with initial 0 value. Choose abritrary positive step.
    } 

    // Choose fit or calculation mode
    printf("\n======================================\n");	      
    if(iparm[38]==0)
    { // Fit
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
            if(param[36]<0.0001) {i=19;} // When chisqr changes less than 1/1000 stop at next iteration. 
            write_io(filename,&COMM,iparm, param,fex,fth,w,fsym);        // Generate output files for next run
            chisqr0=chisqr1;
            write_out(filename,&COMM,iparm,dfact,param,fex,fth,w,fsym,sigmap,dfdc,lookup); // Generate output
            free_vector((double *)fth[0]);
            free_ivector((int *)fsym[0]);
        }
        
    }

    if(iparm[38]>0)
    { // Calculate lines, number given by iparm[38]
		FILE *fp;
        if ((fp = fopen ("freqout_calc.dat", "a")) == NULL) ruserror ("cannot open output file");
		fprintf(fp,"%s",COMM.comment[0]);  // Write title line    
        
		calclines (lookup,param,iparm,dfact,fth,fsym);
        printf("==============================\n");
        print_vector((double *)fth[0],iparm[38]);
		
		if(iparm[39])
        {
            gen_dfdc(dfdc,lookup,param,iparm,dfact,fth,fsym); 
            for(j=1;j<iparm[38]+1;j++)
            {
                printf("%3i ",j); 
				fprintf(fp,"%3i %8.6f",j,((double *)fth[0])[j]); 
                sum=0;
                for(i=1;i<iparm[36]%10+1;i++)
                {
                    sum+= ((double *)dfdc[0])[(i-1)*iparm[38]+j];  
					fprintf(fp,"%5.2f ",((double *)dfdc[0])[(i-1)*iparm[38]+j]);
                    printf("%5.2f ",((double *)dfdc[0])[(i-1)*iparm[38]+j]);
                }
                fprintf(fp,"%5.2f\n",sum);
                printf("%5.2f\n",sum);
            }
            free_vector((double *)dfdc[0]);
        }
        else
		{
			for(j=1;j<iparm[38]+1;j++)
            {
				fprintf(fp,"%3i %8.6f\n",j,((double *)fth[0])[j]); 
			}
		}
	fclose(fp);
    free_vector((double *)fth[0]);
    free_ivector((int *)fsym[0]);        
    }

        if(iparm[38]==-1) // Monte carlo simulation.
    { // Calculate lines
        FILE *fp;
        double chi2;
        int ndata=0;
        
        for(i=1;i<iparm[44]+1;i++)
        {
            if(  ((double *)w[0])[i]!=0.) ndata++; 
        }

        for(i=1;i<=iparm[36]%100;i++)
        {
            mclc[i]=param[i]*(1.-iparm[i]/1000.);
            mcdc[i]=param[i]*2.*iparm[i]/1000.;
        }
         
        time_t t;
        srand((unsigned) time(&t));

        if ((fp = fopen ("chi2out.dat", "a")) == NULL) ruserror ("cannot open output file");
        fprintf(fp,"%s",COMM.comment[0]);  // Write title line
		fprintf(fp,"#Step  chi2     cxx\n");
        fprintf(fp,"#============================================\n");
        fclose(fp);

        for(j=0;j<1000;j=j+1)
        {
            for(i=1;i<=iparm[36]%10;i++)
            {
            param[i]=mclc[i]+mcdc[i]*(double)rand()/(double)RAND_MAX;
            }
        //printf("\n# %i\n",j);
        calclines (lookup,param,iparm,dfact,fth,fsym);
        //print_vector((double *)fth[0],iparm[44]);
        //print_vector((double *)w[0],iparm[44]);
        chi2=fchisq(fex, fth, w,fth,iparm[44],ndata);
        free_vector((double *)fth[0]);
        free_ivector((int *)fsym[0]);
        if ((fp = fopen ("chi2out.dat", "a")) == NULL) ruserror ("cannot open output file");
        printf("Chi2 = %10.6f || ",chi2);
        fprintf(fp,"%4i %10.6f ",j,chi2);
        for(i=1;i<=iparm[36]%100;i++)
            {
            printf("%8.6f ",param[i]);
            fprintf(fp,"%8.6f ",param[i]);
            }
        printf("\n");
        fprintf(fp,"\n");
        //print_vector(param,6);
        fclose(fp);
        }
    }

    
    if(iparm[38]==-2)
    { 
       FILE *fp,*seqmat; 
       double sum; 
       // Calculate lines number given by iparm[18], for a systematic variation of dimensions, at constant density.
      // The dimension to change is given by iparm[37]. 
       b_mass=param[40];
       //b_dens=param[39];
       // Control dimension is 30-32 or 50-52; set in cd.
       cd=0;
       iparm[38]=iparm[44];
       if(iparm[37]>0&&iparm[37]<4) cd=iparm[37]+29;
       if(iparm[37]>3&&iparm[37]<7) cd=iparm[37]+46;
       if(iparm[37]>6&&iparm[37]<10) cd=iparm[37]+26;
       if(cd==0) ruserror("Specify which dimension (1-6) or angle (phi=7; th =8; psi = 9) to vary.\n");
       b_param=param[cd]; // Buffer initial parameter.
       if(b_param==0) ruserror("Variation of dimension or angle can't be done for a parameter equal to zero.\n");
       if ((fp = fopen ("sequence_out.dat", "w")) == NULL) ruserror ("cannot open output file");
       if ((seqmat = fopen ("sequence_matrix.dat", "w")) == NULL) ruserror ("cannot open output file");
       fprintf(seqmat,"%s",COMM.comment[0]);  // Write title line
	   fprintf(fp,"%s",COMM.comment[0]);  // Write title line
	   fprintf(seqmat,"#Parameter|");
       for(i=1;i<=iparm[44];i++)  // Loop on number of lines to calculate.
               {
                fprintf(seqmat,"Line %3i|",i);            
               }
       fprintf(seqmat,"\n");
       for(i=1;i<21;i++)
           {
            param[cd]=b_param/10.*i; // 20 steps from 10% to 200% of the initial value. For angles: advise is to use 90 deg, which calculates 10-200 deg.
            // Line below changes the mass if a dimension is changed. For Euler angles, no change.
            if(iparm[37]<7) {
                            param[40]=b_mass/10.*i;  // For most shapes, mass is linear in the dimensions.
                            switch(iparm[45])
                {
                 // Use mass = density * volume   
                 case 4: param[40] = param[39]/1000.*(param[30]+param[50])*(param[31]+param[51])*(param[32]+param[52])/(6/PI); break;
                 case 8: param[40] = param[39]/1000.*(param[32]*PI/4*(param[30]*param[30]-param[31]*param[31])); break; // Hollow cylinder. OD**2 - ID**2.
                }
                           }
                           
            if(param[40]>0) {
                     calclines (lookup,param,iparm,dfact,fth,fsym);
            printf("=====================================\n");
            printf("= Dimension: %f, Mass: %f\n",param[cd],param[40]);
            printf("=====================================\n\n");
            
            print_vector((double *)fth[0],iparm[44]);
            if(iparm[39]) { gen_dfdc(dfdc,lookup,param,iparm,dfact,fth,fsym); } // Get derivative of frequencies vs moduli.
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
               if(iparm[39]) {free_vector((double *)dfdc[0]);}
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
	}
    
    if(iparm[38]==-3) // Check out error propagation
    { // Calculate lines
	    FILE *fp;
        if ((fp = fopen ("errorout.dat", "a")) == NULL) ruserror ("cannot open output file");
		fprintf(fp,"%s",COMM.comment[0]);  // Write title line
		for(j=1;j<=iparm[36]%100;j++)
        {	
			fprintf(fp,"\n# Deviation of elastic constant %i\n",j);
            printf("\n# Deviation of elastic constant %i\n",j);
            cbuff=param[j];
            for(i=-5;i<=5;i++)
            {
                param[j]=cbuff*(1.+i*1./10000.);
                calclines (lookup,param,iparm,dfact,fth,fsym);
                chisqr0=fchisq(fex, fth, w,fth,iparm[44],iparm[44]);
				fprintf(fp,"%8.6f %8.6f %2i/10000\n",param[j],chisqr0,i);
                printf("%8.6f %8.6f %2i/10000\n",param[j],chisqr0,i);
            }
            param[j]=cbuff;
        }
        printf("\n");
        free_vector((double *)fth[0]);
        free_ivector((int *)fsym[0]);
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
			mclc[i]=param[i]*(1.-iparm[i]/1000.);
			mcdc[i]=param[i]*2.*iparm[i]/1000.;
        }
		if (ngrid==0) ruserror ("no parameter to construct grid on");
        if ((fp = fopen ("chi2out_grid.dat", "a")) == NULL) ruserror ("cannot open output file");
// Handle loops: which cxx to use for grid
// Format output according to ngrid
// Case ngrid=0 do nothing. ngrid directs how many nested loops to do.
		fprintf(fp,"%s",COMM.comment[0]);  // Write title line
		fprintf(fp,"#Step  chi2      c[%i]",iparm[61]);
		if (ngrid>1) fprintf(fp,"      c[%i]",iparm[62]);
		if (ngrid>2) fprintf(fp,"      c[%i]",iparm[63]);
		fprintf(fp,"\n");
        fprintf(fp,"#==============================================================\n");
        fclose(fp);

		i=0;gs=0;
		imax=20*(iparm[61]!=0)+1;
		jmax=20*(iparm[62]!=0)+1;
		kmax=20*(iparm[63]!=0)+1;
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
