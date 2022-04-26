#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include <signal.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#define FL_REAL 0
#define FL_COMPLEX 1

void print_matrix( gsl_matrix_complex * , int , int );
void partial_trace_A(gsl_matrix_complex * , gsl_matrix_complex * , int);
void partial_trace_B(gsl_matrix_complex * , gsl_matrix_complex * , int);
void kronecker_product(gsl_matrix_complex * , gsl_matrix_complex * , gsl_matrix_complex * , int , int , int , int );
void generate_rho0(gsl_matrix_complex * , int , int , int , char *);
void get_random_rho2(gsl_matrix_complex *[7] , gsl_matrix_complex *[4] , int , int , int , int );
void get_random_rhodxd(gsl_matrix_complex * , gsl_matrix_complex *[4], int , int);
double get_trace(gsl_matrix_complex * , int );
//double preselection_criterion(gsl_matrix_complex *[7], int , double *);
//double calculate_p(gsl_matrix_complex *[7] , int , double *);
//double calculate_d(gsl_matrix_complex *[7] , double , double , int , double *, gsl_matrix_complex *[6]);

const gsl_rng_type * T1, * T2;
gsl_rng * rndm, * rndm2;
int NEED_TO_EXIT=0;

void sigterm_handler(int sig){
	if(sig==SIGTERM)
		NEED_TO_EXIT=1;
}


int main(int argc, char * argv[]){
	
	/* Register the signal handler */
	if(signal(SIGTERM,sigterm_handler)==SIG_ERR)
		printf("\nCouldn't catch SIGTERM. Exiting gracelessly.\n");
	
	/* User inputs the number of parties and dimensions as the first and the second argument.*/
	if(argc!=5){
		printf("Wrong/No arguements provided. #n #dim #file #corrections \n");
		return 0;
	}
	
	int dim, np, ncorr;		 /*Dimension and number of parties*/
	int i=0, j=0, k=0, pi=0; /*Counters*/
	double d=99999, prev_d=99999;  /*HSN distance, current and previous*/
	double p; 		 /*Optimum p value for the current run*/
	double fl_pre=-1.;	 /*Preselection flag*/
	int dimn; /*dimension of the n-party matrix*/
	char fname1[26]="pptm\0", fname2[26]="pptm\0"; /* Output files. 1: {n_correction, d}; 2: {{n_correction}, {rho_1}} */
	strcat(fname1, argv[3]);
	strcat(fname2, argv[3]);
	strcat(fname1, "2-n.txt\0");
	strcat(fname2, "2-nr.txt\0");
	gsl_matrix_complex *paramdnxdn[7]; /*0,1,2 store rho0, rho1, rho2. Rest for calculations in the functions*/
	gsl_matrix_complex *paramdxd[4];   /*dxd matrices for calculations in the matrices*/
	gsl_matrix_complex *matdiff[6];   /*Store the differences of the combination of rho0-2*/
	double tr_diff[6];		   /*and in the array, their traces*/
					   /*In order: rho0^2, rho1^2, rho2^2, rho0.rho1, rho0.rho2, rho1.rho2  */
	
	/* Scanning command line arguements */
	sscanf(argv[1], "%d", &np);
	sscanf(argv[2], "%d", &dim);
	sscanf(argv[4], "%d", &ncorr);
	k= strlen(argv[3]);
	dimn = pow(dim, np);
	
	/* Appropriate names for the output files */
	fname1[4+k] = '0'+dim;
	fname2[4+k] = '0'+dim;
	fname1[6+k] = '0'+np;
	fname2[6+k] = '0'+np;

	/*Global Random Number Generator initialization*/
	gsl_rng_env_setup();
	T1 = gsl_rng_ranlxd1;
	T2 = gsl_rng_ranlxd1;
	rndm = gsl_rng_alloc(T1);
	rndm2 = gsl_rng_alloc(T2);
	gsl_rng_set(rndm2, 1234);
	
	/*Memory allocation for everything*/
	for(i=0; i<7; i++){
		paramdnxdn[i]  = gsl_matrix_complex_calloc (dimn,dimn);
	}
	for(i=0; i<6; i++){
		matdiff[i]  = gsl_matrix_complex_calloc (dimn,dimn);
	}
	paramdxd[0] = gsl_matrix_complex_calloc (dim,dim);
	paramdxd[1] = gsl_matrix_complex_calloc (dim,dim);
	paramdxd[2] = gsl_matrix_complex_calloc (dim,dim);	
	paramdxd[3] = gsl_matrix_complex_calloc (dim,dim);
	
	/*Generate or read rho0 from the file*/
	generate_rho0(paramdnxdn[0], dim, dimn, np, argv[3]);
	//printf("\n");

	/***** Generate a dimn identity matrix as rho1******/
/*	gsl_matrix_complex_set_identity(paramdnxdn[1]);*/
/*	gsl_matrix_complex_scale(paramdnxdn[1], gsl_complex_rect(1./(double)dimn, 0.));*/
	
	/*** Or Use the diagonal of rho0 as rho1(better convergence)***/
	gsl_vector_complex_const_view diag_rho0 = gsl_matrix_complex_const_diagonal(paramdnxdn[0]);
	for(i=0; i<dimn; i++)
		gsl_matrix_complex_set(paramdnxdn[1], i,i, gsl_vector_complex_get(&diag_rho0.vector, i));
	/*Update matdiff[0]=rho0-rho1*/
	gsl_matrix_complex_memcpy(matdiff[0],paramdnxdn[0]);
	gsl_matrix_complex_sub(matdiff[0], paramdnxdn[1]);
	
	gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, matdiff[0],matdiff[0], GSL_COMPLEX_ZERO, paramdnxdn[4]);/*Calculate (rho0-rho1)^2*/
	prev_d = get_trace(paramdnxdn[4], dimn); /* and then its trace as the first distance*/
	
	
	/**************************************/
	
/*	If rho1's to be random complex mixed state:*/
/*	get_random_rho2(rho1, dim, FL_COMPLEX);*/
/*	get_random_rho2(rho2, dim, FL_COMPLEX);*/
/*	gsl_matrix_complex_add(rho1, rho2);*/
/*	gsl_matrix_complex_scale(rho1, gsl_complex_rect(0.5, 0.));*/

	//check if the matrix is being read right.
	//print_matrix(paramdnxdn[1], dimn, dimn); printf("\n");

	double darr[ncorr+10];
	
	FILE * o1= NULL;
	o1 = fopen(fname1,"a");
	if(o1==NULL){
		printf("Could not open file.\n");
		return 0;
	}
	//setvbuf(o1, NULL, _IOLBF, 24);
	
	int flag = FL_COMPLEX; //real/complex flag
	int count=0, count_real=0, count_rej=0,interval=0;
	for(i=pi; i<ncorr && d>0.0000000001 && NEED_TO_EXIT==0; i++){
		
		/* Get a rho2 that satisfies the preselection_criterion */
		fl_pre = 0;
		while(fl_pre<=0 || p>=1 || p<0){
			count_rej++;
			get_random_rho2(paramdnxdn, paramdxd, dim, dimn, np, flag);
			/*calculate (rho2-rho1)*/
			gsl_matrix_complex_memcpy(matdiff[1],paramdnxdn[2]);
			gsl_matrix_complex_sub(matdiff[1],paramdnxdn[1]);
			/*calculate (rho2-rho1)^2 and store it's trace (later used to calculate 'p')*/
			gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, matdiff[1],matdiff[1], GSL_COMPLEX_ZERO, matdiff[2]);
			tr_diff[2]=get_trace(matdiff[2],dimn);
			/*calculate (rho0-rho1)(rho2-rho1) and it's trace which is the value of the preselection criterion*/
			gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, matdiff[0],matdiff[1], GSL_COMPLEX_ZERO, paramdnxdn[4]);
			fl_pre=get_trace(paramdnxdn[4], dimn);
			/*from the value of preselection criterion calculate 'p'*/
			p = 1- (fl_pre / tr_diff[2] );
		}
			
		/*calculate rho0+p(rho2-rho1)-rho2, and it's square*/
		/*trace of which is the new distance*/

		gsl_matrix_complex_memcpy(matdiff[3], paramdnxdn[0]);
		gsl_matrix_complex_sub(matdiff[3],paramdnxdn[2]);
		gsl_matrix_complex_scale(matdiff[1],gsl_complex_rect(p,0.));
		gsl_matrix_complex_add(matdiff[3],matdiff[1]);

		gsl_blas_zhemm(CblasLeft, CblasUpper, GSL_COMPLEX_ONE, matdiff[3],matdiff[3], GSL_COMPLEX_ZERO, matdiff[4]);
		d=get_trace(matdiff[4],dimn);

		
		//printf("prev_d=%f, d=%f\n", prev_d, d);
		/* If the new 'd' is less than the 'prev_d', update rho1 */
		/* Otherwise find another rho2, that might reduce distance. */
		if(d >= prev_d) {
			count++;
			if(count>10000){
				flag = (flag + 1) % 2;

				count=0;
			}
			i--;
			continue;
		}else{  /* d < prev_d */
			/* Update rho1 = p*rho1 + (1-p)*rho2 = -matdiff[3]+rho0 */
			gsl_matrix_complex_scale(matdiff[3], gsl_complex_rect( -1   ,0.));
			gsl_matrix_complex_add(matdiff[3],paramdnxdn[0]);
			gsl_matrix_complex_memcpy(paramdnxdn[1],matdiff[3]);
			/*update rho0-rho1*/
			gsl_matrix_complex_memcpy(matdiff[0],paramdnxdn[0]);
			gsl_matrix_complex_sub  (matdiff[0], paramdnxdn[1]);


			count = 0;
			
			if(flag == FL_REAL) {
				count_real++;
				if(count_real>1000){
					flag = FL_COMPLEX;
					count_real=0;
				}
			}
			prev_d = d;
		}
		/* Execution reaches here only if d<prev_d. 					*/
		/* Here we print the values of the rho1 and the distances for each correction.  */
		/* The output is dumped every #interval iterations for the distances, */
		/* and rho1 (Gives best CSS found yet)	                */
		//printf("\n{i=%d, d=%lf} rej:%d\n", i+1, d, count_rej);
		count_rej=0;
		//Output to the files:
		darr[i] = d;
		interval=100;
		if((i+1)%interval==0){
			for(j=i-interval+1; j<i+1; j++){
				fprintf(o1, "{%d, %f}\n",j+1,darr[j]);
			}
			fflush(o1);
		}
		if((i+1)%200==0){
			FILE * o2 = fopen(fname2,"a");
			fprintf(o2, "{{%d},\n{",i+1);
			for(j=0; j<dimn;j++){
				gsl_complex z = gsl_matrix_complex_get(paramdnxdn[1], j, 0);			
				fprintf(o2, "{%lf + I (%lf)", GSL_REAL(z), GSL_IMAG(z));
				for(k=1; k<dimn; k++){
					z = gsl_matrix_complex_get(paramdnxdn[1], j, k);
					fprintf(o2, ", %lf + I (%lf)", GSL_REAL(z), GSL_IMAG(z));
				}
				if(j==dimn-1) fprintf(o2, "}");
				else	  fprintf(o2, "},\n");
			}
			fprintf(o2, "}}\n");
			fflush(o2);
			fclose(o2);	
		}
	}
	/* Memory free */
	if(NEED_TO_EXIT==1)
		printf("\n Flushing and exiting. \n");

	
	if(d<=0.0000000001 && i<2000){
		for(j=i+1-(i+1)%interval; j<i+1; j++){
			fprintf(o1,"{%d, %f}\n",j+1,darr[j]);
		}
		FILE * o2 = fopen(fname2,"a");
		fprintf(o2, "{{%d},\n{",i+1);
		for(j=0; j<dimn;j++){
			gsl_complex z = gsl_matrix_complex_get(paramdnxdn[1], j, 0);			
			fprintf(o2, "{%lf + I (%lf)", GSL_REAL(z), GSL_IMAG(z));
			for(k=1; k<dimn; k++){
				z = gsl_matrix_complex_get(paramdnxdn[1], j, k);
				fprintf(o2, ", %lf + I (%lf)", GSL_REAL(z), GSL_IMAG(z));
			}
			if(j==dimn-1) fprintf(o2, "}");
			else	  fprintf(o2, "},\n");
		}
		fprintf(o2, "}}\n");
		fflush(o2);
		fclose(o2);	
	}
	fflush(o1);
	fclose(o1);
	
	for(i=0; i<7; i++){
		gsl_matrix_complex_free(paramdnxdn[i]);
	}
	for(i=0; i<6;i++){
		gsl_matrix_complex_free(matdiff[i]);
	}
	gsl_matrix_complex_free(paramdxd[0]);
	gsl_matrix_complex_free(paramdxd[1]);
	gsl_matrix_complex_free(paramdxd[2]);
	gsl_matrix_complex_free(paramdxd[3]);

	gsl_rng_free(rndm);
	gsl_rng_free(rndm2);

	return 0;
}

void generate_rho0(gsl_matrix_complex * rho0, int dim, int dimn, int np , char * nm){ //Maximally entangled n qudit state
	char name[20] = "ppt\0";
	strcat(name, nm);
	strcat(name,".dat\0");
	FILE * lp = fopen(name,"r");
	gsl_matrix_complex_fscanf(lp, rho0);
	fclose(lp);
	return;
}

void print_matrix( gsl_matrix_complex * mat, int m, int n){
	int i=0, j=0;
	gsl_complex z;
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			z = gsl_matrix_complex_get (mat, i, j);
			printf("%lf + I (%lf) ", GSL_REAL(z), GSL_IMAG(z));
		}
		printf("\n");
	}
	return;
}

void kronecker_product (gsl_matrix_complex * A, gsl_matrix_complex * B, gsl_matrix_complex * R, int nA, int mA, int nB, int mB){
/* na <= nb && ma <= mb */
	int i=0, j=0;
	for(i=0;i<nA;i++){
		for(j=0;j<mA;j++){
			gsl_complex aij = gsl_matrix_complex_get(A, i,j);
			if( GSL_IMAG(aij)!=0 || GSL_REAL(aij)!=0 ){
				gsl_matrix_complex_view sR = gsl_matrix_complex_submatrix(R, i*nB,j*mB, nB,mB);
				gsl_matrix_complex_memcpy(&sR.matrix, B);
				gsl_matrix_complex_scale(&sR.matrix, aij);
			}
		}
	}
	return;	
}

void get_random_rho2(gsl_matrix_complex * p1[7], gsl_matrix_complex * p2[4], int dim, int dimn, int np, int flag){
	gsl_matrix_complex * rhoA = p2[0];
	gsl_matrix_complex * intr = p1[3];
	gsl_matrix_complex * rho2 = p1[2];
	gsl_matrix_complex_set_zero(p1[2]);
	gsl_matrix_complex_view rhoB = gsl_matrix_complex_submatrix(intr, 0,0, dim,dim);
	gsl_matrix_complex_view res;
	get_random_rhodxd(rhoA, p2, dim, flag);
	get_random_rhodxd(&rhoB.matrix, p2, dim, flag);
	int d[np], i;
	
	d[0] = dim;
	for(i=1; i<np; i++)
		d[i] = d[i-1]*dim;
	res = gsl_matrix_complex_submatrix(rho2, 0,0, d[1],d[1]);
	gsl_matrix_complex_set_zero(&res.matrix);
	kronecker_product(rhoA, &rhoB.matrix, &res.matrix, dim, dim, d[0], d[0]);
	
	for(i=1; i<np-1; i++){
		rhoB = gsl_matrix_complex_submatrix(intr, 0,0, d[i],d[i]);
		gsl_matrix_complex_memcpy(&rhoB.matrix, &res.matrix);
		res = gsl_matrix_complex_submatrix(rho2, 0,0, d[i+1], d[i+1]);
		gsl_matrix_complex_set_zero(&res.matrix);
		get_random_rhodxd(rhoA, p2, dim, flag);
		kronecker_product(rhoA, &rhoB.matrix, &res.matrix, dim, dim, d[i], d[i]);
	}
	return;
}

void get_random_rhodxd(gsl_matrix_complex * rho4, gsl_matrix_complex * p2[4], int dim, int flag){
	gsl_matrix_complex * a = p2[2];
	int i, pr1=0, pr2=0; 
	double rn1r=0., rn1i=0., rn2=0., rn22=0.;
	for (i=0;i<dim;i++){
		rn2  = gsl_rng_uniform(rndm2);
		rn22 = gsl_rng_uniform(rndm2);
		rn1r = gsl_rng_uniform( rndm);
		rn1i = gsl_rng_uniform( rndm);
		if(rn2>0.5) pr1=1; else pr1=-1;
		if(rn22>0.5) pr2=1; else pr2=-1;		
		gsl_matrix_complex_set(a, i,0, gsl_complex_rect(pr1*sqrt((-2.0)*log(rn1r)), pr2*flag*sqrt((-2.0)*log(rn1i))) );
	}
	gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, a, a, GSL_COMPLEX_ZERO, rho4);
	gsl_matrix_complex_scale( rho4, gsl_complex_rect(1/get_trace(rho4, dim), 0.) );
	return;
}

double get_trace(gsl_matrix_complex * rho, int n){
	int i;
	gsl_complex tr = GSL_COMPLEX_ZERO;
	for(i=0; i<n; i++){
		tr = gsl_complex_add(tr,gsl_matrix_complex_get(rho,i,i));
	}
	return GSL_REAL(tr);
}

/**/
/**/
/**/
/*END*/
