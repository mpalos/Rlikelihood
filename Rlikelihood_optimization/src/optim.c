/*
 * TODO Include verbose levels
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
# include <unistd.h>

#include <gsl/gsl_vector.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_multifit_nlin.h>

#include <time.h>
#include <sys/timeb.h>
#include <sys/types.h>

#include "lnLfunctions.h"

/******************************************************************************************/
/* FUNCTION PROTOTYPES */
/******************************************************************************************/
void printResult(char model, int iter,
		gsl_vector * xVal, double lik, double grdSizeAcc,
		int status);
/******************************************************************************************/

#define handle_error(msg) \
	do {perror(msg); exit(EXIT_FAILURE);} while (0)

// Global variables:
int nObs = 0;
int *PSobs = NULL;
int *SSobs = NULL;
int *DSobs = NULL;
int *PNobs = NULL;
int *SNobs = NULL;
int *DNobs = NULL;
//double *x_init = NULL;
double step_size, tol, epsabs;
double maxSame, maxIter;
//double Ls = 3000, Ln = 7000, n = 24;
double cubTol, cubZero, cubSInf;
int cubMaxEval;
double *x_init = NULL;
double *ss_init = NULL;
//long nData = 0;
//int startPoint = -1;

// Global variables:
double *x_params = NULL;
int nParam = -1;
int nRerun = 0;

void printResult(char model, int iter,
		gsl_vector * xVal, double lik, double grdSizeAcc,
		int status){

	  printf("%d;",iter);
	  printf("%.20f;",exp(gsl_vector_get(xVal, 0)));
	  printf("%.20f;",exp(gsl_vector_get(xVal, 1)));
	  printf("%.20f;",exp(gsl_vector_get(xVal, 2)));
	  printf("%.20f;",gsl_vector_get(xVal, 3));
	  printf("%.20f;",gsl_vector_get(xVal, 4));
	  if(model == 'B' || model == 'C'){
		  printf("%.20f;",1 / (1 + exp(-1*gsl_vector_get (xVal, 5))));
		  printf("%.20f;",exp(gsl_vector_get (xVal, 6)));
	  }
	  printf("%.20f;",exp(gsl_vector_get(xVal, 7)));
	  printf("%.20f;",exp(gsl_vector_get(xVal, 8)));
	  printf("%.20f;", lik);
	  printf("%.20f;",grdSizeAcc);
	  printf("%d\n",status);

}

double lnL_f(const gsl_vector * x, void *params) {

	double lnL = 0.0;
	int j;
	const char model = ((struct params_s *) params)->model;

	double theta = gsl_vector_get(x, 0);
	double lambda = gsl_vector_get(x, 1);
	double r = gsl_vector_get(x, 2);
	double smean = gsl_vector_get(x, 3);
	double smax = gsl_vector_get(x, 4);
	double Pb = gsl_vector_get(x, 5);
	double Sb = gsl_vector_get(x, 6);
	double shapeBeta = gsl_vector_get(x, 7);
	double shapeMu = gsl_vector_get(x, 8);

//	  /**** DEBUG *****/
//	  printf("x: %f %f %f %f %f %f %f %f %f \n"
//			  ,gsl_vector_get(x,0),gsl_vector_get(x,1),gsl_vector_get(x,2),gsl_vector_get(x,3)
//			  ,gsl_vector_get(x,4),gsl_vector_get(x,5),gsl_vector_get(x,6),gsl_vector_get(x,7)
//			  ,gsl_vector_get(x,8));
//	  /****************/

	// PARAMETERS
	if(model == 'A'){
		theta 			= exp(gsl_vector_get(x, 0));
		lambda 			= exp(gsl_vector_get(x, 1));
		r 				= exp(gsl_vector_get(x, 2));
		smean 			= gsl_vector_get(x, 3);
		smax 			= gsl_vector_get(x, 4);
		Pb				= 0;
		Sb				= 0;
		shapeBeta 		= exp(gsl_vector_get(x, 7));
		shapeMu 		= exp(gsl_vector_get(x, 8));
	}

	else if(model == 'B' || model == 'C'){
		theta 			= exp(gsl_vector_get(x, 0));
		lambda 			= exp(gsl_vector_get(x, 1));
		r 				= exp(gsl_vector_get(x, 2));
		smean 			= gsl_vector_get(x, 3);
		smax		= -cubZero;
		Pb			= 1 / (1 + exp(-1*gsl_vector_get (x, 5)));
		Sb			= exp(gsl_vector_get (x, 6));
		shapeBeta 		= exp(gsl_vector_get(x, 7));
		shapeMu 		= exp(gsl_vector_get(x, 8));
	}

//	/**** DEBUG ******/
//	printf("x: %f %f %f %f %f %f %f %f %f \n",theta,lambda,r,smean,smax,Pb,Sb,shapeBeta,shapeMu);
//	/*****************/

	for (j = 0; j < ((struct params_s *) params)->nObs; j++) {
		lnL += lnL_MKS_PhiS_PhiMu(model,
				((struct params_s *) params)->Ln,
				((struct params_s *) params)->Ls,
				((struct params_s *) params)->n,
				((struct params_s *) params)->PSobs[j],
				((struct params_s *) params)->PNobs[j],
				((struct params_s *) params)->SNobs[j],
				((struct params_s *) params)->SSobs[j],
				((struct params_s *) params)->DNobs[j],
				((struct params_s *) params)->DSobs[j],
				theta, lambda, r, smean, smax, Pb, Sb, shapeBeta, shapeMu,
				cubTol, cubMaxEval, cubZero, cubSInf);

		/****** DEBUG *****/
//		printf("obs: %d %d %d %d %d %d %f\n",
//				((struct params_s *) params)->PSobs[j],
//				((struct params_s *) params)->PNobs[j],
//				((struct params_s *) params)->SNobs[j],
//				((struct params_s *) params)->SSobs[j],
//				((struct params_s *) params)->DNobs[j],
//				((struct params_s *) params)->DSobs[j],
//				lnL_MKS_PhiS_PhiMu(model,
//								((struct params_s *) params)->Ln,
//								((struct params_s *) params)->Ls,
//								((struct params_s *) params)->n,
//								((struct params_s *) params)->PSobs[j],
//								((struct params_s *) params)->PNobs[j],
//								((struct params_s *) params)->SNobs[j],
//								((struct params_s *) params)->SSobs[j],
//								((struct params_s *) params)->DNobs[j],
//								((struct params_s *) params)->DSobs[j],
//								theta, lambda, r, smean, smax, Pb, Sb, shapeBeta, shapeMu,
//								cubTol, cubMaxEval, cubZero, cubSInf));
//
		/*****************/

	}

	const double y = lnL * (-1);

	return y;
}

int lnL_f1(const gsl_vector * x, void *params, gsl_vector * f) {

	gsl_vector_set(f, 0, lnL_f(x, params));

	return GSL_SUCCESS;

}

int lnL_df(const gsl_vector *x, void *params, gsl_vector *df) {

	//HERE IS THE DEFERENTIAL USING FINITE DIFFERENCES (multifit_fdfsolver)
	gsl_multifit_function_fdf fdf;
	fdf.f = &lnL_f1;
	fdf.df = NULL;
	fdf.fdf = NULL;
	fdf.n = 1;
	fdf.p = 9;
	fdf.params = params;

	const char model = ((struct params_s *) params)->model;

	gsl_matrix *J = gsl_matrix_alloc(1, 9);
	J->size1 = 1;
	J->size2 = 9;

	gsl_vector *f = gsl_vector_alloc(1);

	gsl_vector_set(f, 0, lnL_f(x, params));

	int status = gsl_multifit_fdfsolver_dif_df(x, &fdf, f, J);

	if (status) {
		printf("ERROR #%d", status);
		return status;
	}

	gsl_vector_set(df, 0, gsl_matrix_get(J, 0, 0));
	gsl_vector_set(df, 1, gsl_matrix_get(J, 0, 1));
	gsl_vector_set(df, 2, gsl_matrix_get(J, 0, 2));
	gsl_vector_set(df, 3, gsl_matrix_get(J, 0, 3));
	gsl_vector_set(df, 4, gsl_matrix_get(J, 0, 4));
	if(model == 'B' || model == 'C'){
		gsl_vector_set(df, 5, gsl_matrix_get(J, 0, 5));
		gsl_vector_set(df, 6, gsl_matrix_get(J, 0, 6));
	}
	gsl_vector_set(df, 7, gsl_matrix_get(J, 0, 7));
	gsl_vector_set(df, 8, gsl_matrix_get(J, 0, 8));

	return GSL_SUCCESS;

}

int lnL_fdf(const gsl_vector * x, void *params, double * f, gsl_vector *df) {

	*f = lnL_f(x, params);
	lnL_df(x, params, df);
	return GSL_SUCCESS;
}

int optim_lnL_derivatives(const gsl_multimin_fdfminimizer_type *T, void *params,
		gsl_vector *x, gsl_vector *res, double *f, double *grad, int *outIter) {

	unsigned int iter = 0;
	unsigned int countSame = 0;
	int status;
	gsl_multimin_fdfminimizer *s;

	gsl_vector *prev = gsl_vector_alloc(9);

	gsl_multimin_function_fdf my_func;
	my_func.n = 9;
	my_func.f = lnL_f;
	my_func.df = lnL_df;
	my_func.fdf = lnL_fdf;
	my_func.params = params;

	const char model = ((struct params_s *) params)->model;

	s = gsl_multimin_fdfminimizer_alloc(T, 9);

	gsl_multimin_fdfminimizer_set(s, &my_func, x, step_size, tol); //ORIGINAL

//	  /**** DEBUG *****/
//	  printf("x: %f %f %f %f %f %f %f %f %f \n"
//			  ,gsl_vector_get(x,0),gsl_vector_get(x,1),gsl_vector_get(x,2),gsl_vector_get(x,3)
//			  ,gsl_vector_get(x,4),gsl_vector_get(x,5),gsl_vector_get(x,6),gsl_vector_get(x,7)
//			  ,gsl_vector_get(x,8));
//
//	  /***************/

	//PRINT HEADER IF NECESSARY
//	printf ("startPoint;nData;iter;"
//			"theta;lambda;r;smean;smax;shapeBeta;shapeMu;"
//			"f;grad;status\n"
//	 );

//	  /**** DEBUG *****/
//	  printf("x: %f %f %f %f %f %f %f %f %f \n"
//			  ,gsl_vector_get(x,0),gsl_vector_get(x,1),gsl_vector_get(x,2),gsl_vector_get(x,3)
//			  ,gsl_vector_get(x,4),gsl_vector_get(x,5),gsl_vector_get(x,6),gsl_vector_get(x,7)
//			  ,gsl_vector_get(x,8));
//	  /****************/

	double mod_grad = 0.0;
	int n_dim = 2;

	//PRINT START POINT
	mod_grad = pow(gsl_vector_get(s->gradient, 0), n_dim)
			+ pow(gsl_vector_get(s->gradient, 1), n_dim)
			+ pow(gsl_vector_get(s->gradient, 2), n_dim)
			+ pow(gsl_vector_get(s->gradient, 3), n_dim)
			+ pow(gsl_vector_get(s->gradient, 4), n_dim)
			+ pow(gsl_vector_get(s->gradient, 7), n_dim)
			+ pow(gsl_vector_get(s->gradient, 8), n_dim);

	  if(model == 'B' || model == 'C'){
		  mod_grad += pow(gsl_vector_get(s->gradient, 5), n_dim)
					+ pow(gsl_vector_get(s->gradient, 6), n_dim);
	  }

	  printResult(model,iter,s->x,s->f,pow(mod_grad, 1.0 / n_dim),-2);

	gsl_vector_memcpy(prev, x);

	do {
		iter++;
		status = gsl_multimin_fdfminimizer_iterate(s);
//		  /**** DEBUG *****/
//		  printf("x: %f %f %f %f %f %f %f %f %f \n"
//				  ,gsl_vector_get(x,0),gsl_vector_get(x,1),gsl_vector_get(x,2),gsl_vector_get(x,3)
//				  ,gsl_vector_get(x,4),gsl_vector_get(x,5),gsl_vector_get(x,6),gsl_vector_get(x,7)
//				  ,gsl_vector_get(x,8));
//		  /****************/

		if (status) {
			// IF IS NECESSARY TO PRINT THE ERRORS ...
			//if (status==27){
			//	printf("ERROR #%d: iteration is not making progress towards solution\n",status);
			//}

			printResult(model,iter,s->x,s->f,pow(mod_grad, 1.0 / n_dim),status);

			break;
		}

		status = gsl_multimin_test_gradient(s->gradient, epsabs);
		// IF IS NECESSARY TO PRINT THE MIN ...
		//if (status == GSL_SUCCESS){
		//	printf ("Minimum found at:\n");
		//}

		mod_grad = pow(gsl_vector_get(s->gradient, 0), n_dim)
				+ pow(gsl_vector_get(s->gradient, 1), n_dim)
				+ pow(gsl_vector_get(s->gradient, 2), n_dim)
				+ pow(gsl_vector_get(s->gradient, 3), n_dim)
				+ pow(gsl_vector_get(s->gradient, 4), n_dim)
				+ pow(gsl_vector_get(s->gradient, 7), n_dim)
				+ pow(gsl_vector_get(s->gradient, 8), n_dim);

		  if(model == 'B' || model == 'C'){
			  mod_grad += pow(gsl_vector_get(s->gradient, 5), n_dim)
						+ pow(gsl_vector_get(s->gradient, 6), n_dim);
		  }

		  printResult(model,iter,s->x,s->f,pow(mod_grad, 1.0 / n_dim),status);

		  //TEST IF THE ITERATION STUCKS ON SAME ANSWER AFTER maxSame TIMES
		  if(gsl_vector_equal(prev,s->x)){ countSame++; }
		  else { countSame =0;}

		  if(countSame >= maxSame){
			  printResult(model,iter,s->x,s->f,pow(mod_grad, 1.0 / n_dim),27);
			  break;
			}

		gsl_vector_memcpy(prev, s->x);

	} while (status == GSL_CONTINUE && iter < maxIter);

	//VERIFY IF AFTER maxIter ITERATIONS IT COULDN'T FIND SOLUTION
	if (iter >= maxIter) {
		printResult(model,iter,s->x,s->f,pow(mod_grad, 1.0 / n_dim),27);
	}

	gsl_vector_memcpy(res, s->x);
	f[0] = s->f;
	grad[0] = pow(mod_grad, 1.0 / n_dim);
	outIter[0] = iter;

	gsl_multimin_fdfminimizer_free(s);


	return status;

}

int optim_lnL_simplex(const gsl_multimin_fminimizer_type *T, void *params, gsl_vector *x, gsl_vector *ss,
		gsl_vector *res, double *f, double *outSize, int *outIter) {

	  unsigned int n = 9;

	  gsl_multimin_fminimizer *s = NULL;

	  unsigned int iter = 0;
	  int status;
	  double size = 0.0;
	  const char model = ((struct params_s *) params)->model;

	  gsl_multimin_function minex_func = {&lnL_f, n, params};

//	  /**** DEBUG *****/
//	  printf("x: %f %f %f %f %f %f %f %f %f \n"
//			  ,gsl_vector_get(x,0),gsl_vector_get(x,1),gsl_vector_get(x,2),gsl_vector_get(x,3)
//			  ,gsl_vector_get(x,4),gsl_vector_get(x,5),gsl_vector_get(x,6),gsl_vector_get(x,7)
//			  ,gsl_vector_get(x,8));
//
//	  printf("ss: %f %f %f %f %f %f %f %f %f \n"
//			  ,gsl_vector_get(ss,0),gsl_vector_get(ss,1),gsl_vector_get(ss,2),gsl_vector_get(ss,3)
//			  ,gsl_vector_get(ss,4),gsl_vector_get(ss,5),gsl_vector_get(ss,6),gsl_vector_get(ss,7)
//			  ,gsl_vector_get(ss,8));
//
//	  /***************/

	  //Use x_params variable
	  gsl_vector_set(ss, 0, log(1 + (gsl_vector_get(ss,0) / exp(gsl_vector_get(x,0)))));
	  gsl_vector_set(ss, 1, log(1 + (gsl_vector_get(ss,1) / exp(gsl_vector_get(x,1)))));
	  gsl_vector_set(ss, 2, log(1 + (gsl_vector_get(ss,2)   / exp(gsl_vector_get(x,2)))));
	  gsl_vector_set(ss, 3, gsl_vector_get(ss,3));
	  gsl_vector_set(ss, 4, gsl_vector_get(ss,4));

	  if (model == 'B' || model == 'C'){
		  gsl_vector_set(x,5,1 / (1 + exp(-1*gsl_vector_get (x, 5)))); // Back to the Pb original value
	  	  gsl_vector_set(ss, 5, log((gsl_vector_get(x,5) + gsl_vector_get(ss,5))/(1-(gsl_vector_get(x,5) + gsl_vector_get(ss,5)))) - log(gsl_vector_get(x,5)/(1-gsl_vector_get(x,5))));
	  	  gsl_vector_set(x,5,log(gsl_vector_get(x,5)/(1-gsl_vector_get(x,5)))); // Now the Pb transformation

	  	  gsl_vector_set(ss, 6, gsl_vector_get(ss,6));
	  }

	  gsl_vector_set(ss, 7, log(1 + (gsl_vector_get(ss,7)  / exp(gsl_vector_get(x,7)))));
	  gsl_vector_set(ss, 8, log(1 + (gsl_vector_get(ss,8)  / exp(gsl_vector_get(x,8)))));


//	  /**** DEBUG *****/
//	  printf("ss: %f %f %f %f %f %f %f %f %f \n"
//			  ,gsl_vector_get(ss,0),gsl_vector_get(ss,1),gsl_vector_get(ss,2),gsl_vector_get(ss,3)
//			  ,gsl_vector_get(ss,4),gsl_vector_get(ss,5),gsl_vector_get(ss,6),gsl_vector_get(ss,7)
//			  ,gsl_vector_get(ss,8));
//
//	  /***************/

	  gsl_set_error_handler_off();

	  s = gsl_multimin_fminimizer_alloc (T, n);
	  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	  //PRINT START POINT
	  size = gsl_multimin_fminimizer_size (s);
//	  /**** DEBUG *****/
//	  printf("size: %f\n",s->size);
//	  /****************/


	  printResult(model,iter,s->x,lnL_f(x,params),size,-2);

	  	unsigned int var = 0;
	  	unsigned int iterPerRun = 0;
		for (var = 0; var <= nRerun; var++) {
		  	iterPerRun = 0; //set iterPerRun
			do
			  {
				iter++;
				iterPerRun++;
				status = gsl_multimin_fminimizer_iterate(s);

				if (status){ //print iteration for error status
					if(var!=nRerun){
						status = status + (100 * (nRerun - var)); // To indicate intermediate execs
					}
					printResult(model,iter,s->x,s->fval,size,status);

					if(var!=nRerun){ //Reset rerun
						gsl_multimin_fminimizer_set (s, &minex_func, s->x, ss); //set the simplex again!!!
						//iterPerRun = 0; //set iterPerRun
					}

					break;
				}

				size = gsl_multimin_fminimizer_size (s);
				status = gsl_multimin_test_size (size, epsabs);

				if(status == 0 && var!=nRerun){ //Force to status if do not reach the # of re-executions
					status = status + (100 * (nRerun - var)); // To indicate intermediate execs
					gsl_multimin_fminimizer_set (s, &minex_func, s->x, ss); //set the simplex again!!!
					//iterPerRun = 0; //set iterPerRun
				}

				printResult(model,iter,s->x,s->fval,size,status);


			  }
			while (status == GSL_CONTINUE && iterPerRun < maxIter);

			//PRINT STATUS IF IT HAVEN'T FOUND THE SOLUTION AFTER maxIter ITERATIONS PER RUN
			if(iterPerRun >= maxIter){
				status = 27 + (100 * (nRerun - var));
				printResult(model,iter,s->x,s->fval,size,status);
			}
		}

	gsl_vector_memcpy(res, s->x);
	f[0] = s->fval;
	outSize[0] = size;
	outIter[0] = iter;

	//gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	  return status;
}

int runOptimization_der(char model, char method
		,int nObs
		,int *obsPs, int *obsPn, int *obsSn, int *obsSs, int *obsDn, int *obsDs
		,double *startPoint
		,double *constants
		,double *methodParams
		,double *cubParams
		,double *outRes
		,double *outf
		,double *outGrad
		,int *outIter){

	double Ls, Ln, n;
	Ls = constants[0];
	Ln = constants[1];
	n = constants[2];

	epsabs = methodParams[0];
	maxIter = methodParams[1];
	step_size = methodParams[2];
	tol = methodParams[3];
	maxSame = methodParams[4];

//	/***** DEBUG *******/
//	printf("params: %f %f %f %f %f\n",methodParams[0],methodParams[1],methodParams[2],methodParams[3],methodParams[4]);
//	printf("params: %f %f %f %f %f\n",epsabs,maxIter,step_size,tol,maxSame);
//
//	/*******************/

	cubZero = cubParams[0];
	cubSInf = cubParams[1];
	cubMaxEval = cubParams[2];
	cubTol = cubParams[3];

	PSobs = obsPs;
	SSobs = obsSs;
	DSobs = obsDs;
	PNobs = obsPn;
	SNobs = obsSn;
	DNobs = obsDn;

	params par;
	par.Ls = Ls;
	par.Ln = Ln;
	par.n = n;
	par.nObs = nObs;
	par.PSobs = PSobs;
	par.SSobs = SSobs;
	par.DSobs = DSobs;
	par.PNobs = PNobs;
	par.SNobs = SNobs;
	par.DNobs = DNobs;
	par.method = method;
	par.model = model;

	//startpoint memory allocation
	x_init = malloc(sizeof(double) * 9); if (x_init == NULL) {handle_error("malloc: x_init");}
	// x_init <- startPoint
	if(model == 'A'){
		x_init[0] = startPoint[0];		// theta
		x_init[1] = startPoint[1];		// lambda
		x_init[2] = startPoint[2];		// r
		x_init[3] = startPoint[3];		// smean
		x_init[4] = startPoint[4];		// smax
		x_init[7] = startPoint[7];		// shapeBeta
		x_init[8] = startPoint[8];		// shapeMu
	}

	else if(model == 'B' || model == 'C'){
		x_init[0] = startPoint[0]; 		// theta
		x_init[1] = startPoint[1];		// lambda
		x_init[2] = startPoint[2];		// r
		x_init[3] = startPoint[3];		// smean
		x_init[4] = startPoint[4];		// smax
		x_init[5] = startPoint[5];		// Pb
		x_init[6] = startPoint[6];   	// Sb
		x_init[7] = startPoint[7];		// shapeBeta
		x_init[8] = startPoint[8];		// shapeMu
	}


	gsl_vector_view x = gsl_vector_view_array(x_init, 9);

	// Define optimization algorithm/ run the optimizer

	const gsl_multimin_fdfminimizer_type *Tfdfmin;

	switch (method) {
		case 'B':
			Tfdfmin = gsl_multimin_fdfminimizer_vector_bfgs2;
			break;
		case 'P':
			Tfdfmin = gsl_multimin_fdfminimizer_conjugate_pr;
			break;
		case 'F':
			Tfdfmin = gsl_multimin_fdfminimizer_conjugate_fr;
			break;
		case 'D':
			Tfdfmin = gsl_multimin_fdfminimizer_steepest_descent;
			break;
		default:
			fprintf(stderr,"Error: method %c not valid",method);
			return EXIT_FAILURE;
			break;
	}


	int status;

	gsl_vector *res = gsl_vector_alloc(9);

	status = optim_lnL_derivatives(Tfdfmin, &par, &x.vector, res, outf, outGrad, outIter);

//	/************ DEBUG *********/
//		printf("res: %f %f %f %f %f %f %f %f %f %f %f %f\n"
//				,gsl_vector_get(res,0)
//				,gsl_vector_get(res,1)
//				,gsl_vector_get(res,2)
//				,gsl_vector_get(res,3)
//				,gsl_vector_get(res,4)
//				,gsl_vector_get(res,5)
//				,gsl_vector_get(res,6)
//				,gsl_vector_get(res,7)
//				,gsl_vector_get(res,8)
//				,gsl_vector_get(addRes,0)
//				,gsl_vector_get(addRes,1)
//				,gsl_vector_get(addRes,2));
//
//	/************       *********/

	//get results
	outRes[0] = exp(gsl_vector_get(res,0));
	outRes[1] = exp(gsl_vector_get(res,1));
	outRes[2] = exp(gsl_vector_get(res,2));
	outRes[3] = gsl_vector_get(res,3);
	outRes[4] = gsl_vector_get(res,4);
	if(model == 'A'){
		outRes[5] = 0;
		outRes[6] = 0;
	}
	else if(model == 'B' || model == 'C'){
		outRes[5] = 1 / (1 + exp(-1*gsl_vector_get (res, 5)));
		outRes[6] = exp(gsl_vector_get(res,6));
	}
	outRes[7] = exp(gsl_vector_get(res,7));
	outRes[8] = exp(gsl_vector_get(res,8));

	return status;


}

int runOptimization_wder(char model, char method
		,int nObs
		,int *obsPs, int *obsPn, int *obsSn, int *obsSs, int *obsDn, int *obsDs
		,double *startPoint
		,double *constants
		,double *stepSize
		,double *methodParams
		,double *cubParams
		,double *outRes
		,double *outf
		,double *outSize
		,int *outIter){


	double Ls, Ln, n;
	unsigned int nT = 9;

	Ls = constants[0];
	Ln = constants[1];
	n = constants[2];

	epsabs = methodParams[0];
	maxIter = methodParams[1];
	nRerun = methodParams[2];

	cubZero = cubParams[0];
	cubSInf = cubParams[1];
	cubMaxEval = cubParams[2];
	cubTol = cubParams[3];

	PSobs = obsPs;
	SSobs = obsSs;
	DSobs = obsDs;
	PNobs = obsPn;
	SNobs = obsSn;
	DNobs = obsDn;

	params par;
	par.Ls = Ls;
	par.Ln = Ln;
	par.n = n;
	par.nObs = nObs;
	par.PSobs = PSobs;
	par.SSobs = SSobs;
	par.DSobs = DSobs;
	par.PNobs = PNobs;
	par.SNobs = SNobs;
	par.DNobs = DNobs;
	par.method = method;
	par.model = model;

	//startpoint memory allocation
	x_init = malloc(sizeof(double) * nT); if (x_init == NULL) {handle_error("malloc: x_init");}
	// x_init <- startPoint
	//if(model == 'A'){
		x_init[0] = startPoint[0];	// theta
		x_init[1] = startPoint[1];	// lambda
		x_init[2] = startPoint[2];	// r
		x_init[3] = startPoint[3];	// smean
		x_init[4] = startPoint[4];	// smax
		x_init[5] = startPoint[5];	// Pb
		x_init[6] = startPoint[6];	// Sb
		x_init[7] = startPoint[7];	// shapeBeta
		x_init[8] = startPoint[8];	// shapeMu
	//}

//	else if(model == 'B' || model == 'C'){
//		x_init[0] = startPoint[0]; 	// theta
//		x_init[1] = startPoint[1];	// lambda
//		x_init[2] = startPoint[2];	// r
//		x_init[3] = startPoint[3];	// smean
//		x_init[4] = startPoint[4];	// smax
//		x_init[5] = startPoint[5];	// Pb
//		x_init[6] = startPoint[6];   	// Sb
//		x_init[7] = startPoint[7];	// shapeBeta
//		x_init[8] = startPoint[8];	// shapeMu
//	}

	//stepSize memory allocation
	ss_init = malloc(sizeof(double) * nT); if (x_init == NULL) {handle_error("malloc: ss_init");}
	// x_init <- startPoint
//	if(model == 'A'){
		ss_init[0] = stepSize[0];		// theta
		ss_init[1] = stepSize[1];		// lambda
		ss_init[2] = stepSize[2];		// r
		ss_init[3] = stepSize[3];		// smean
		ss_init[4] = stepSize[4];		// smax
		ss_init[5] = stepSize[5];		// Pb
		ss_init[6] = stepSize[6];		// Sb
		ss_init[7] = stepSize[7];		// shapeBeta
		ss_init[8] = stepSize[8];		// shapeMu
//	}

//	else if(model == 'B' || model == 'C'){
//		ss_init[0] = stepSize[0]; 	// theta
//		ss_init[1] = stepSize[1];		// lambda
//		ss_init[2] = stepSize[2];		// r
//		ss_init[3] = stepSize[3];		// smean
//		ss_init[4] = stepSize[4];		// smax
//		ss_init[5] = stepSize[5];		// Pb
//		ss_init[6] = stepSize[6];   	// Sb
//		ss_init[7] = stepSize[7];		// shapeBeta
//		ss_init[8] = stepSize[8];		// shapeMu
//	}

	gsl_vector_view x = gsl_vector_view_array(x_init, nT);
	gsl_vector_view ss = gsl_vector_view_array(ss_init, nT);

	/****** DEBUG ********/
	// Print x_init
//	printf("%f %f %f %f %f %f %f %f %f\n"
//			,x_init[0],x_init[1],x_init[2],x_init[3],x_init[4],x_init[5],x_init[6],x_init[7],x_init[8]);
//
//	printf("%f %f %f %f %f %f %f %f %f\n"
//			,stepSize[0],stepSize[1],stepSize[2],stepSize[3],stepSize[4],stepSize[5],stepSize[6],stepSize[7],stepSize[8]);
	/*********************/

	// Define optimization algorithm/ run the optimizer

	const gsl_multimin_fminimizer_type *Tfmin;

	if (method == 'S'){
		Tfmin = gsl_multimin_fminimizer_nmsimplex2;
	}
	else {
		fprintf(stderr,"Error: method %c not valid",method);
		return EXIT_FAILURE;
	}

	int status;

	gsl_vector *res = gsl_vector_alloc(9);

	status = optim_lnL_simplex(Tfmin, &par, &x.vector, &ss.vector, res, outf, outSize, outIter);

	//get results
	outRes[0] = exp(gsl_vector_get(res,0));
	outRes[1] = exp(gsl_vector_get(res,1));
	outRes[2] = exp(gsl_vector_get(res,2));
	outRes[3] = gsl_vector_get(res,3);
	outRes[4] = gsl_vector_get(res,4);
	if(model == 'A'){
		outRes[5] = 0;
		outRes[6] = 0;
	}
	else if(model == 'B' || model == 'C'){
		outRes[5] = 1 / (1 + exp(-1*gsl_vector_get (res, 5)));
		outRes[6] = exp(gsl_vector_get(res,6));
	}
	outRes[7] = exp(gsl_vector_get(res,7));
	outRes[8] = exp(gsl_vector_get(res,8));

	return status;
}
