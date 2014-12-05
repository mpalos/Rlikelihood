/*
 * Rlikelihood.c
 *
 *  Created on: Mar 3, 2014
 *      Author: MAPF
 */

/*
 *
 * R CMD SHLIB *.c -o libRlikelihood.so
 *
 * */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "lnLfunctions.h"

SEXP likelihood(SEXP model, SEXP Ln, SEXP Ls, SEXP n
		,SEXP obsPs, SEXP obsPn, SEXP obsSn, SEXP obsSs, SEXP obsDn, SEXP obsDs
		,SEXP Mutheta, SEXP lambda, SEXP r, SEXP smean, SEXP smax, SEXP Pb, SEXP Sb
		,SEXP shapeBeta, SEXP shapeMu
		,SEXP cubTol, SEXP cubMaxEval, SEXP cubZero, SEXP cubSInf){

	SEXP lnL;

	// GETTING THE FUNCTION PARAMETERS
	int nObs, j; //INTERNAL COUNTER
	char modelID;
	const char *model_;
	double *obsPs_, *obsPn_, *obsSn_, *obsSs_, *obsDn_, *obsDs_;
	double *Ln_, *Ls_, *n_, *Mutheta_, *lambda_, *r_, *smean_, *smax_, *Pb_, *Sb_, *shapeBeta_, *shapeMu_;
	double *cubTol_, *cubMaxEval_, *cubZero_, *cubSInf_;
	double *res;

	//Coerce Vectors to real
	PROTECT(obsPs = coerceVector(obsPs,REALSXP));
	PROTECT(obsPn = coerceVector(obsPn,REALSXP));
	PROTECT(obsSn = coerceVector(obsSn,REALSXP));
	PROTECT(obsSs = coerceVector(obsSs,REALSXP));
	PROTECT(obsDn = coerceVector(obsDn,REALSXP));
	PROTECT(obsDs = coerceVector(obsDs,REALSXP));

	// Now we can attribute to internal C

	obsPs_ = REAL(obsPs);
	obsPn_ = REAL(obsPn);
	obsSn_ = REAL(obsSn);
	obsSs_ = REAL(obsSs);
	obsDn_ = REAL(obsDn);
	obsDs_ = REAL(obsDs);

	int i = 0;
	for(i = 0; i < LENGTH(model); i++){
		model_ = CHAR(STRING_ELT(model, i));
	}

	if(strcmp(model_,"model_A") == 0){ modelID = 'A';	}
	else if(strcmp(model_,"model_B") == 0){	modelID = 'B'; }
	else if(strcmp(model_,"model_C") == 0){	modelID = 'C'; }
	else {error("Invalid model: %s",model_);}

	// Parameters/Constants
	Ln_ = REAL(Ln);
	Ls_ = REAL(Ls);
	n_ = REAL(n);
	Mutheta_ = REAL(Mutheta);
	lambda_ = REAL(lambda);
	r_ = REAL(r);
	smean_ = REAL(smean);
	smax_ = REAL(smax);
	Pb_ = REAL(Pb);
	Sb_ = REAL(Sb);
	shapeBeta_ = REAL(shapeBeta);
	shapeMu_ = REAL(shapeMu);

	//Cubature extra arguments
	cubTol_ = REAL(cubTol);
	cubMaxEval_ = REAL(cubMaxEval);
	cubZero_ = REAL(cubZero);
	cubSInf_ = REAL(cubSInf);

	// WARNING TREATMENT: all Obs must have the same length
	nObs = length(obsPs);
	if(nObs != length(obsPn)) {
		warning("There are some vectors which haven't the same length. Using the shortest vector number of elements");
		if(nObs > length(obsPn)) { nObs = length(obsPn); }
	}
	if(nObs != length(obsSn)) {
		warning("There are some vectors which haven't the same length. Using the shortest vector number of elements");
		if(nObs > length(obsSn)) { nObs = length(obsSn); }
	}
	if(nObs != length(obsSs)) {
		warning("There are some vectors which haven't the same length. Using the shortest vector number of elements");
		if(nObs > length(obsSs)) { nObs = length(obsSs); }
	}
	if(nObs != length(obsDn)) {
		warning("There are some vectors which haven't the same length. Using the shortest vector number of elements");
		if(nObs > length(obsDn)) { nObs = length(obsDn); }
	}
	if(nObs != length(obsDs)) {
		warning("There are some vectors which haven't the same length. Using the shortest vector number of elements");
		if(nObs > length(obsDs)) { nObs = length(obsDs); }
	}

	//ALLOC lnL
	PROTECT(lnL = allocVector(REALSXP,1));
	res = REAL(lnL);
	res[0] = 0.0;
	for(j=0; j< nObs; j++){

		res[0] += lnL_MKS_PhiS_PhiMu(modelID,Ln_[0], Ls_[0], n_[0]
				,obsPs_[j],obsPn_[j],obsSn_[j],obsSs_[j],obsDn_[j],obsDs_[j]
				,Mutheta_[0],lambda_[0],r_[0],smean_[0],smax_[0], Pb_[0], Sb_[0]
				,shapeBeta_[0],shapeMu_[0]
				,cubTol_[0],cubMaxEval_[0],cubZero_[0],cubSInf_[0]);

	}

	UNPROTECT(7);
	return lnL;

}

SEXP calcEPs(SEXP Ls, SEXP n
		,SEXP Mutheta, SEXP r){

	SEXP resEPs;

	// GETTING THE FUNCTION PARAMETERS
	double *Ls_, *Mutheta_, *r_, *n_;
	double *res;

	Ls_ = REAL(Ls);
	n_ = REAL(n);
	Mutheta_ = REAL(Mutheta);
	r_ = REAL(r);

	//ALLOC resEPs
	PROTECT(resEPs = allocVector(REALSXP,1));
	res = REAL(resEPs);
	res[0] = 0.0;

	res[0] += EPs(Ls_[0], r_[0], Mutheta_[0], n_[0]);

	UNPROTECT(1);
	return resEPs;

}

SEXP calcESs (SEXP Ls, SEXP n
				  ,SEXP Mutheta){

	SEXP resESs;

	// GETTING THE FUNCTION PARAMETERS
	double *Ls_, *Mutheta_, *n_;
	double *res;

	Ls_ = REAL(Ls);
	n_ = REAL(n);
	Mutheta_ = REAL(Mutheta);

	//ALLOC resESs
	PROTECT(resESs = allocVector(REALSXP,1));
	res = REAL(resESs);
	res[0] = 0.0;

	res[0] += ESs(Ls_[0], Mutheta_[0], n_[0]);

	UNPROTECT(1);
	return resESs;

}

SEXP calcEDs (SEXP Ls, SEXP n
				  ,SEXP Mutheta, SEXP lambda ){

	SEXP resEDs;

	// GETTING THE FUNCTION PARAMETERS
	double *Ls_, *Mutheta_, *lambda_, *n_;
	double *res;

	Ls_ = REAL(Ls);
	n_ = REAL(n);
	Mutheta_ = REAL(Mutheta);
	lambda_ = REAL(lambda);

	//ALLOC resESs
	PROTECT(resEDs = allocVector(REALSXP,1));
	res = REAL(resEDs);
	res[0] = 0.0;

	res[0] += EDs(Ls_[0], lambda_[0], Mutheta_[0], n_[0]);

	UNPROTECT(1);
	return resEDs;

}

SEXP calcEPn_OverPhiS_adapt(SEXP model, SEXP Ln, SEXP n
		,SEXP Mutheta, SEXP r, SEXP smean, SEXP smax, SEXP Pb, SEXP Sb
		,SEXP shapeBeta
		,SEXP cubTol, SEXP cubMaxEval, SEXP cubZero, SEXP cubSInf){

	SEXP resEPn;

	// GETTING THE FUNCTION PARAMETERS
	const char *model_;
	char modelID;
	double *Ln_, *n_, *Mutheta_, *r_, *smean_, *smax_, *Pb_, *Sb_, *shapeBeta_;
	double *cubTol_, *cubMaxEval_, *cubZero_, *cubSInf_;
	double *res;

	int i = 0;
	for(i = 0; i < LENGTH(model); i++){
		model_ = CHAR(STRING_ELT(model, i));
	}

	if(strcmp(model_,"model_A") == 0){ modelID = 'A';	}
	else if(strcmp(model_,"model_B") == 0){	modelID = 'B'; }
	else if(strcmp(model_,"model_C") == 0){	modelID = 'C'; }
	else {error("Invalid model: %s",model_);}

	// Parameters/Constants
	Ln_ = REAL(Ln);
	n_ = REAL(n);
	Mutheta_ = REAL(Mutheta);
	r_ = REAL(r);
	smean_ = REAL(smean);
	smax_ = REAL(smax);
	Pb_ = REAL(Pb);
	Sb_ = REAL(Sb);
	shapeBeta_ = REAL(shapeBeta);

	//Cubature extra arguments
	cubTol_ = REAL(cubTol);
	cubMaxEval_ = REAL(cubMaxEval);
	cubZero_ = REAL(cubZero);
	cubSInf_ = REAL(cubSInf);

	//ALLOC resESs
	PROTECT(resEPn = allocVector(REALSXP,1));
	res = REAL(resEPn);
	res[0] = 0.0;

	res[0] += calc_EPn_OverPhiS_adapt(modelID,Ln_[0], n_[0]
	                                   , Mutheta_[0], r_[0], smean_[0], smax_[0] ,Pb_[0], Sb_[0]
	                                   , shapeBeta_[0]
	                                   ,cubTol_[0],cubMaxEval_[0],cubZero_[0],cubSInf_[0]);
	UNPROTECT(1);
	return resEPn;


}

SEXP calcESn_OverPhiS_adapt(SEXP model, SEXP Ln, SEXP n
		,SEXP Mutheta, SEXP smean, SEXP smax, SEXP Pb, SEXP Sb
		,SEXP shapeBeta
		,SEXP cubTol, SEXP cubMaxEval, SEXP cubZero, SEXP cubSInf){

	SEXP resESn;

	// GETTING THE FUNCTION PARAMETERS
	char modelID;
	const char *model_;
	double *Ln_, *n_, *Mutheta_, *smean_, *smax_, *Pb_, *Sb_, *shapeBeta_;
	double *cubTol_, *cubMaxEval_, *cubZero_, *cubSInf_;
	double *res;

	int i = 0;
	for(i = 0; i < LENGTH(model); i++){
		model_ = CHAR(STRING_ELT(model, i));
	}

	if(strcmp(model_,"model_A") == 0){ modelID = 'A';	}
	else if(strcmp(model_,"model_B") == 0){	modelID = 'B'; }
	else if(strcmp(model_,"model_C") == 0){	modelID = 'C'; }
	else {error("Invalid model: %s",model_);}

	// Parameters/Constants
	Ln_ = REAL(Ln);
	n_ = REAL(n);
	Mutheta_ = REAL(Mutheta);
	smean_ = REAL(smean);
	smax_ = REAL(smax);
	Pb_ = REAL(Pb);
	Sb_ = REAL(Sb);
	shapeBeta_ = REAL(shapeBeta);

	//Cubature extra arguments
	cubTol_ = REAL(cubTol);
	cubMaxEval_ = REAL(cubMaxEval);
	cubZero_ = REAL(cubZero);
	cubSInf_ = REAL(cubSInf);

	//ALLOC resESs
	PROTECT(resESn = allocVector(REALSXP,1));
	res = REAL(resESn);
	res[0] = 0.0;

	res[0] += calc_ESn_OverPhiS_adapt(modelID, Ln_[0], n_[0]
	                                   , Mutheta_[0], smean_[0], smax_[0], Pb_[0], Sb_[0]
	                                   , shapeBeta_[0]
	                                   ,cubTol_[0],cubMaxEval_[0],cubZero_[0],cubSInf_[0]);
	UNPROTECT(1);
	return resESn;


}

SEXP calcEDn_OverPhiS_adapt(SEXP model, SEXP Ln, SEXP n
		,SEXP Mutheta, SEXP lambda, SEXP r ,SEXP smean, SEXP smax, SEXP Pb, SEXP Sb
		,SEXP shapeBeta
		,SEXP cubTol, SEXP cubMaxEval, SEXP cubZero, SEXP cubSInf){

	SEXP resEDn;

	// GETTING THE FUNCTION PARAMETERS
	const char *model_;
	char modelID;
	double *Ln_, *n_, *Mutheta_, *lambda_, *r_, *smean_, *smax_, *Pb_, *Sb_, *shapeBeta_;
	double *cubTol_, *cubMaxEval_, *cubZero_, *cubSInf_;
	double *res;

	int i = 0;
	for(i = 0; i < LENGTH(model); i++){
		model_ = CHAR(STRING_ELT(model, i));
	}

	if(strcmp(model_,"model_A") == 0){ modelID = 'A';	}
	else if(strcmp(model_,"model_B") == 0){	modelID = 'B'; }
	else if(strcmp(model_,"model_C") == 0){	modelID = 'C'; }
	else {error("Invalid model: %s",model_);}

	// Parameters/Constants
	Ln_ = REAL(Ln);
	n_ = REAL(n);
	lambda_ = REAL(lambda);
	r_ = REAL(r);
	Mutheta_ = REAL(Mutheta);
	smean_ = REAL(smean);
	smax_ = REAL(smax);
	Pb_ = REAL(Pb);
	Sb_ = REAL(Sb);
	shapeBeta_ = REAL(shapeBeta);

	//Cubature extra arguments
	cubTol_ = REAL(cubTol);
	cubMaxEval_ = REAL(cubMaxEval);
	cubZero_ = REAL(cubZero);
	cubSInf_ = REAL(cubSInf);

	//ALLOC resESs
	PROTECT(resEDn = allocVector(REALSXP,1));
	res = REAL(resEDn);
	res[0] = 0.0;

	res[0] += calc_EDn_OverPhiS_adapt(modelID, Ln_[0], n_[0]
	                                  , Mutheta_[0], lambda_[0], r_[0], smean_[0], smax_[0], Pb_[0], Sb_[0]
	                                  , shapeBeta_[0]
	                                  ,cubTol_[0],cubMaxEval_[0],cubZero_[0],cubSInf_[0]);
	UNPROTECT(1);
	return resEDn;


}
