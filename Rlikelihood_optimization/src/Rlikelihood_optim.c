/*
 * Rlikelihood.c
 *
 *  Created on: Mar 3, 2014
 *      Author: MAPF
 */

/* R CMD SHLIB Rlikelihood.c
 *
 * For the cluster execution:
 * MAKEFLAGS="PKG_CFLAGS=-I/home/mfranco/addlibs/include/"  R CMD SHLIB *.c -o libRlikelihood_optimization.so -L/home/mfranco/addlibs/lib -lgslcblas -lgsl
 *
 * */

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

#include "lnLfunctions.h"

SEXP optim_der(SEXP model, SEXP method
		,SEXP obsPs, SEXP obsPn, SEXP obsDs, SEXP obsDn, SEXP obsSs, SEXP obsSn
		,SEXP startPoint
		,SEXP constants
		,SEXP methodParams
		,SEXP cubParams){

	SEXP status;
	SEXP outRes, outf, outGrad, outIter;

	//List names
	SEXP listNames, list;
	char *names[7] = {"model","method","iterations","result","f","grad","status"};

	//C internals
	int nObs; //INTERNAL COUNTER
	char modelID, methodID;
	const char *model_ = NULL, *method_ = NULL;
	int *obsPs_, *obsPn_, *obsSn_, *obsSs_, *obsDn_, *obsDs_;
	double *startPoint_;
	double *constants_;
	double *methodParams_;
	double *cubParams_;
	double *outRes_, *outf_, *outGrad_;
	int *outIter_;
	int *res;

	//Coerce Vectors to int
	PROTECT(obsPs = coerceVector(obsPs,INTSXP));
	PROTECT(obsPn = coerceVector(obsPn,INTSXP));
	PROTECT(obsSn = coerceVector(obsSn,INTSXP));
	PROTECT(obsSs = coerceVector(obsSs,INTSXP));
	PROTECT(obsDn = coerceVector(obsDn,INTSXP));
	PROTECT(obsDs = coerceVector(obsDs,INTSXP));

	PROTECT(constants = coerceVector(constants,REALSXP));
	PROTECT(methodParams = coerceVector(methodParams,REALSXP));
	PROTECT(cubParams = coerceVector(cubParams,REALSXP));
	PROTECT(startPoint = coerceVector(startPoint,REALSXP));

	// Now we can attribute to internal C int

	obsPs_ = INTEGER(obsPs);
	obsPn_ = INTEGER(obsPn);
	obsSn_ = INTEGER(obsSn);
	obsSs_ = INTEGER(obsSs);
	obsDn_ = INTEGER(obsDn);
	obsDs_ = INTEGER(obsDs);

	constants_ = REAL(constants);
	methodParams_ = REAL(methodParams);
	cubParams_ = REAL(cubParams);
	startPoint_ = REAL(startPoint);

	int i = 0;
	for(i = 0; i < LENGTH(model); i++){
		model_ = CHAR(STRING_ELT(model, i));
	}

	if(strcmp(model_,"model_A") == 0){ modelID = 'A';	}
	else if(strcmp(model_,"model_B") == 0){	modelID = 'B'; }
	else if(strcmp(model_,"model_C") == 0){	modelID = 'C'; }
	else {error("Invalid model: %s",model_);}

	for(i = 0; i < LENGTH(method); i++){
		method_ = CHAR(STRING_ELT(method, i));
	}

	if(strcmp(method_,"bfgs") == 0){ methodID = 'B';	}
	else if(strcmp(method_,"conj_pr") == 0){ methodID = 'P'; }
	else if(strcmp(method_,"conj_fr") == 0){ methodID = 'F'; }
	else if(strcmp(method_,"steep_desc") == 0){ methodID = 'D'; }
	else {error("Invalid method: %s",method_);}

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


	PROTECT(outRes = allocVector(REALSXP,9));
	outRes_ = REAL(outRes);

	PROTECT(outf = allocVector(REALSXP,1));
	outf_ = REAL(outf);

	PROTECT(outGrad = allocVector(REALSXP,1));
	outGrad_ = REAL(outGrad);

	PROTECT(outIter = allocVector(INTSXP,1));
	outIter_ = INTEGER(outIter);

	PROTECT(status = allocVector(INTSXP,1));
	res = INTEGER(status);
	res[0] = -1;

	res[0] = runOptimization_der(modelID,methodID
			,nObs,obsPs_,obsPn_,obsSn_,obsSs_,obsDn_,obsDs_
			,startPoint_
			,constants_
			,methodParams_,
			cubParams_
			,outRes_
			,outf_
			,outGrad_
			,outIter_);

//	/****** DEBUG *****/
//	Rprintf("outRes: %f %f %f %f %f %f %f %f %f\n"
//			, outRes_[0], outRes_[1], outRes_[2], outRes_[3], outRes_[4], outRes_[5]
//			, outRes_[6], outRes_[7], outRes_[8]);
//	Rprintf("f: %f", outf_[0]);
//	/******************/

	//CREATE A LIST THAT CONTAINS THE RESULTS
	//LIST NAMES
	PROTECT(listNames = allocVector(STRSXP,7));

	for(i=0 ; i<7 ; i++){
		SET_STRING_ELT(listNames,i,mkChar(names[i]));
	}

	//LIST VALUES
	PROTECT(list = allocVector(VECSXP,7));

	//ATTACHING VALUES INTO LIST
	SET_VECTOR_ELT(list,0,model); 		//model
	SET_VECTOR_ELT(list,1,method); 		//method
	SET_VECTOR_ELT(list,2,outIter); 	//# of iterations
	SET_VECTOR_ELT(list,3,outRes); 		//result vector
	SET_VECTOR_ELT(list,4,outf); 		//function value
	SET_VECTOR_ELT(list,5,outGrad); 	//gradient value
	SET_VECTOR_ELT(list,6,status); 		//status

	setAttrib(list, R_NamesSymbol, listNames); // Include list names

	UNPROTECT(17);

	return list;

}

SEXP optim_wder(SEXP model, SEXP method
		,SEXP obsPs, SEXP obsPn, SEXP obsDs, SEXP obsDn, SEXP obsSs, SEXP obsSn
		,SEXP startPoint
		,SEXP constants
		,SEXP stepSize
		,SEXP methodParams
		,SEXP cubParams){

	SEXP status;
	SEXP outRes, outf, outSize, outIter;

	//List names
	SEXP listNames, list;
	char *names[7] = {"model","method","iterations","result","f","size","status"};

	//C internals
	int nObs; //INTERNAL COUNTER
	char modelID, methodID;
	const char *model_ = NULL, *method_ = NULL;
	int *obsPs_, *obsPn_, *obsSn_, *obsSs_, *obsDn_, *obsDs_;
	double *startPoint_;
	double *constants_;
	double *stepSize_;
	double *methodParams_;
	double *cubParams_;
	double *outRes_, *outf_, *outSize_;
	int *outIter_;
	int *res;

	//Coerce Vectors to int
	PROTECT(obsPs = coerceVector(obsPs,INTSXP));
	PROTECT(obsPn = coerceVector(obsPn,INTSXP));
	PROTECT(obsSn = coerceVector(obsSn,INTSXP));
	PROTECT(obsSs = coerceVector(obsSs,INTSXP));
	PROTECT(obsDn = coerceVector(obsDn,INTSXP));
	PROTECT(obsDs = coerceVector(obsDs,INTSXP));

	PROTECT(constants = coerceVector(constants,REALSXP));
	PROTECT(methodParams = coerceVector(methodParams,REALSXP));
	PROTECT(cubParams = coerceVector(cubParams,REALSXP));
	PROTECT(startPoint = coerceVector(startPoint,REALSXP));
	PROTECT(stepSize = coerceVector(stepSize,REALSXP));


	// Now we can attribute to internal C int

	obsPs_ = INTEGER(obsPs);
	obsPn_ = INTEGER(obsPn);
	obsSn_ = INTEGER(obsSn);
	obsSs_ = INTEGER(obsSs);
	obsDn_ = INTEGER(obsDn);
	obsDs_ = INTEGER(obsDs);

	constants_ = REAL(constants);
	methodParams_ = REAL(methodParams);
	cubParams_ = REAL(cubParams);
	startPoint_ = REAL(startPoint);
	stepSize_ = REAL(stepSize);

	int i = 0;
	for(i = 0; i < LENGTH(model); i++){
		model_ = CHAR(STRING_ELT(model, i));
	}

	if(strcmp(model_,"model_A") == 0){ modelID = 'A';	}
	else if(strcmp(model_,"model_B") == 0){	modelID = 'B'; }
	else if(strcmp(model_,"model_C") == 0){	modelID = 'C'; }
	else {error("Invalid model: %s",model_);}

	for(i = 0; i < LENGTH(method); i++){
		method_ = CHAR(STRING_ELT(method, i));
	}

	if(strcmp(method_,"simplex") == 0){ methodID = 'S';	}
	else {error("Invalid method: %s",method_);}

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


//	/****** DEBUG ********/
//	// Print  startPoint_ and stepSize_
//	printf("startPoint_: %f %f %f %f %f %f %f %f %f\n"
//			,startPoint_[0],startPoint_[1],startPoint_[2],startPoint_[3],startPoint_[4],startPoint_[5],startPoint_[6],startPoint_[7],startPoint_[8]);
//
//	printf("%f %f %f %f %f %f %f %f %f\n"
//			,stepSize_[0],stepSize_[1],stepSize_[2],stepSize_[3],stepSize_[4],stepSize_[5],stepSize_[6],stepSize_[7],stepSize_[8]);
//	/*********************/

	PROTECT(outRes = allocVector(REALSXP,9));
	outRes_ = REAL(outRes);

	PROTECT(outf = allocVector(REALSXP,1));
	outf_ = REAL(outf);

	PROTECT(outSize = allocVector(REALSXP,1));
	outSize_ = REAL(outSize);

	PROTECT(outIter = allocVector(INTSXP,1));
	outIter_ = INTEGER(outIter);

	PROTECT(status = allocVector(INTSXP,1));
	res = INTEGER(status);
	res[0] = -1;

	res[0] = runOptimization_wder(modelID,methodID
			,nObs,obsPs_,obsPn_,obsSn_,obsSs_,obsDn_,obsDs_
			,startPoint_
			,constants_
			,stepSize_
			,methodParams_,
			cubParams_
			,outRes_
			,outf_
			,outSize_
			,outIter_);

	//CREATE A LIST THAT CONTAINS THE RESULTS
	//LIST NAMES
	PROTECT(listNames = allocVector(STRSXP,7));

	for(i=0 ; i<7 ; i++){
		SET_STRING_ELT(listNames,i,mkChar(names[i]));
	}

	//LIST VALUES
	PROTECT(list = allocVector(VECSXP,7));

	//ATTACHING VALUES INTO LIST
	SET_VECTOR_ELT(list,0,model); 		//model
	SET_VECTOR_ELT(list,1,method); 		//method
	SET_VECTOR_ELT(list,2,outIter); 	//# of iterations
	SET_VECTOR_ELT(list,3,outRes); 		//result vector
	SET_VECTOR_ELT(list,4,outf); 		//function value
	SET_VECTOR_ELT(list,5,outSize); 	//gradient value
	SET_VECTOR_ELT(list,6,status); 		//status

	setAttrib(list, R_NamesSymbol, listNames); // Include list names

	UNPROTECT(18);

	return list;

}
