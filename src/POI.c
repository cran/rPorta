/*
 * This file is part of RPorta. For licensing and copyright information
 * please see the file COPYING in the root directory of this
 * distribution or contact <robin.nunkesser@tu-dortmund.de>.
 */ 

#include <R.h>
#include <Rinternals.h>  
#include <Rdefines.h>
//#include <R_ext/Rdynload.h>

#include "POI.h"	
#include "arith.h"
#include "porta.h"
//#include "valid.h"
#include "common.h"
#include "inout.h"
#include "limits.h"


//#define _DEBUG_TODO
//#define _DEBUG
//#define DEBUG_R_RUECKGABE
//#ifdef COMPILE_AS_LIB
//#endif
	// ******************************************************************************
	// ***********************************************************************Eingabe
 
//int dimMatrix  [2];// =  { 5, 3};


 
int precision = 1000;


//double * matrix;
//double ** versagenPunkte  = 0;
	// ******************************************************************************	




	 
long loRound(double x) {
	if (x > 0) 
		return (long) x + 0.5; 
	return (long) x - 0.5;
} 




void classifyPoints(RAT * matrix, int * dimMatrix, int * failureInd, RAT ** outFailurePoints, RAT ** outSuccessPoints, int *failureCount, int *successCount){
	/*  @param:  *matrix ^= pointer to double array which contains the points to be examinated (int row major) -> (k1,...,kn),(k1,...,kn)...
	*            *dimMatrix ^= pointer to int array which contains dimension of Matrix 
                               note: (dimMatrix[0] = Count of points, dimMatrix[1] = dimension of points
				 *failureInd ^= pointer to an int array (length = count of points) failurePoints[i] = 1 means that point i in Matrix is
									  a failure point
				 *outFailurePoints ^= pointer which will be setted to an new array that contains all FailurePoints
				 *outSuccessPoints ^= pointer which will be setted to an new array that contains all Erfolgspunkte
				 *failureCount ^= stores the number of failurePoints
	 */
	int i;
	(*failureCount)=0;
	for (i = 0; i < dimMatrix[0]; i++){
		if (failureInd[i]){
			(*failureCount) += 1;
		}
	}
	*successCount = (dimMatrix[0]-(*failureCount));
	(*outFailurePoints) = (RAT *) malloc ((*failureCount) * dimMatrix[1] * sizeof(RAT));
	(*outSuccessPoints) = (RAT *) malloc ((*successCount) * dimMatrix[1] * sizeof(RAT));
	
	int fC,eC,j;
	fC = eC = 0;
	
	for (i = 0; i < dimMatrix[0]; i++) {
		if (failureInd[i]){
			for (j = 0; j < dimMatrix[1]; j++){
				(*outFailurePoints)[fC * dimMatrix[1] + j].num = matrix[i * dimMatrix[1] + j].num;
				(*outFailurePoints)[fC * dimMatrix[1] + j].den.i = matrix[i * dimMatrix[1] + j].den.i;
			}
			fC++;
		} else {
			for (j = 0; j < dimMatrix[1]; j++){
				(*outSuccessPoints)[eC * dimMatrix[1] + j].num = matrix[i * dimMatrix[1] + j].num;
				(*outSuccessPoints)[eC * dimMatrix[1] + j].den.i = matrix[i * dimMatrix[1] + j].den.i;
			}
			eC++;
		}
	}
}
		
	
	

void I_RAT_Init(RAT * a, long pnum, long pi){
	a->num = pnum;
	a->den.i = pi;
}

void I_RAT_Kuerzen(RAT * a){
	int _gcd = igcd(a->num,a->den.i);
	a->num /= _gcd;
	a->den.i /= _gcd;
}

int I_RAT_Vergleich(RAT * a , RAT * b){
	
	int tempA = a->num * b->den.i;
	int tempB = b->num * a->den.i;
	if (tempA < tempB) return -1;
	if (tempA > tempB) return  1;
	return  0;
}

//###############################################################################
// ##########################################################################TEST

void constructPolyCones(RAT * matSucc, int * dimSucc, RAT * matFail, int * dimFail,
						RAT ** polyCones){
	/*  constructPolyCones constructs Poly Cones represented by an array of Vectors
	 *  the Vectors will be stored in an 3 dimensional Array of Rational Numbers (RAT) 
	 *  the First Dimension selects a PolyCone, the Secound Dimension the Vector of the 
	 *  PolyCone and the Third Dimension indicates the Component of the Vector of the PolyCone
	 *  thus polyCones[2][1][3] indicates the forth component of the second Vector of the third 
	 *  polyCone (note: the indication starts with zero)
	 *
	 *  for each Failure Point one PolyCone will be constructed. each PolyCone has 
	 *  dimSucc[0] vectors which will be constructed by subtraction the given FailurePoint 
	 *  of each SucsessPoint. 
	 *  hence we get dimFail[0] polyCones with each dimSucc[0] vectors of length 
	 *  dimSucc[0] = dimFail[0]
	 *
	 *  @matSucc a pointer to an array which contains the Success Points
	 *  @dimSucc a pointer to an array of length 2 which contains the number and the Dimension
	 *			 of the Success Points
	 *  @matFail a pointer to an array which contains the Failure Points
	 */
	
	
	// polyCones wird in dieser Funktion initialisiert
	//polyCones = 0;
	//TODO Warum dieses = 0? Das führt zu einem Absturz
	(*polyCones) = (RAT *) RATallo(*polyCones,0,dimSucc[0] * dimSucc[1] * dimFail[0]);
	
	int vecNrV =   dimFail[0];      
	int vecNrR =   dimSucc[0];
	int dimVec =   dimSucc[1];
	
	
	//RAT  a;
	int c, v, d ;
	for (c = 0;  c < vecNrV; c++){
		for (v = 0; v < vecNrR; v++){
			for (d = 0; d < dimVec; d++){
				/*a.num = (int)((matFail[c*dimVec +d]-matSucc[v*dimVec +d ]) * 1000000000 + .5);                        a.den.i = 1000000000;                         	                         	                       
				 (*polyCones)[ c * (vecNrR*dimVec) +  v * dimVec + d ] = a; */
				//RAT fail, succ;
				
				//I_RAT_Init(&fail,loRound(matFail[c*dimVec +d] * precision),precision);
				//I_RAT_Init(&succ,loRound(matSucc[v*dimVec +d] * precision),precision);
			
				RAT a;
				I_RAT_sub(matFail[c*dimVec +d],matSucc[v*dimVec +d],&a);
				
				(*polyCones)[ c * (vecNrR*dimVec) +  v * dimVec + d ] = a;

				
			}
		}
	}
	
	/* hier steht:
	 *      -- *1000 000 000 um Rundungsvehler 9 Nachkomastele zu senken;
	 *      -- +0.5 damit Runden richtig funktioniert:  Bp: (int)0.7 = 0 was wir nicht wollen;  
	 *         aber (int)(0.7 +0.5) = 1;                                  
	 */  
}	

int testPointInHyperPlane(RAT * halfspace, RAT * point, int pointDimension){
	/*
	 @param		*plane ^= pointer to an inequality represented by an array of RAT's 
	 (k1 , k2, .... , k_pointDimension, b)	 which can be read as 
	 k1 * x1 + k2 * x2 + ... + k_pointDimension * x_pointDimension <= b
					*point    ^= pointer to an Array of RAT's which represent the point to test 
	 (must have an Dimension of pointDimension)
	 @function:  checks whether the point satisfies the inequality  
	 */
//---------------------------------------------------------------------------------
	//printf("Enter inHyper?");	
	#ifdef _DEBUG 
		printf("Ueberpruefe Gleichung:\n");
		printIneq(halfspace,pointDimension,1);
		printf("Mit Punkt:\n");
		printar(point,pointDimension,0,1);
		//if (within) printf("Punkt erfuellt die Gleichung\n\n");
		//else printf("Punkt erfuellt nicht die Gleichung\n\n");
	#endif
	 
	int i;
	RAT leftSide;
	leftSide.num = 0; leftSide.den.i = 1;
	for (i = 0; i < pointDimension; i++){
		RAT addition;
		I_RAT_mul(*(halfspace+i),*(point+i),&addition);
		I_RAT_add(addition,leftSide,&leftSide);
	}
	int within = I_RAT_Vergleich(&leftSide,halfspace+pointDimension) != 1;
//---------------------------------------------------------------------------------	
//	#ifdef _DEBUG 
//		printf("Ueberpruefe Gleichung:\n");
//		printIneq(halfspace,pointDimension,1);
//		printf("Mit Punkt:\n");
//		printar(point,pointDimension,0,1);
//		if (within) printf("Punkt erfuellt die Gleichung\n\n");
//		else printf("Punkt erfuellt nicht die Gleichung\n\n");
//	#endif
	return within;
	
} 


int testPointInPolyCone(RAT * polyCone, RAT * point, int faceCount, int pointDimension){	
	/*
	 @param:  *polyCone ^= pointer to an array of Inequalities stored as follows: 
	 IEQ1 , IEQ2 , ..., IEQ_faceCount
	 an Inequality is represented as an Array of pointDimension + 1
	 RAT's (k1 , k2, .... , k_pointDimension, b)  which can be read as 
	 k1 * x1 + k2 * x2 + ... + k_pointDimension * x_pointDimension <= b
	 faceCount ^= number of Inequalities in polyCone 
	 *point    ^= pointer to an Array of RAT's which represent the point to test 
	 (must have an Dimension of pointDimension)
					pointDimension ^= dimension of the point 
	 @function:  check whether the point satisfies all inequalities
	 */
	int i,pointIn;
	pointIn = 1; // nothing is tested yet -> pointIn is true 
	RAT * aktInequality = polyCone;
	for (i = 0; (i <faceCount) && (pointIn); i++){
		pointIn = pointIn && testPointInHyperPlane(aktInequality,point,pointDimension);
		aktInequality=aktInequality+(pointDimension+1);
	}
	return pointIn;
}  




int testPointInUnionOfPolyCones(RAT ** polyCones, RAT * coneTranslation, int * polyConesFacesCount, int polyCount, RAT * point, int pointDimension, int * ieq){
	int i,j;
	#ifdef _DEBUG
	printf("------------- Ueberprfe Punkt------------:\n");
	printar(point,pointDimension,0,1);
	printf("------------- Translationspunkte (versagenspunkte) ------------:\n");
	printar(coneTranslation,pointDimension,0,polyCount);
	#endif
	
	for (i = 0; i < polyCount;i++){
		RAT tempTrans[pointDimension];
		for (j = 0; j < pointDimension; j++){ /*a-b*/
			I_RAT_sub(point[j],coneTranslation[i*pointDimension+j],&tempTrans[j]);
		}
		if (testPointInPolyCone(polyCones[i],(RAT*) (&tempTrans),polyConesFacesCount[i],pointDimension)){
			
			//'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
			//printf("SET index-------: %i\n",(*ieq));			
			(*ieq) =  i;
			//printf("DONE index-------: %i\n",(*ieq));	
					
			//''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
			
			#ifdef _DEBUG
			printf("------------- Ende der Ueberprfung von Punkt ------------:\n");
			printar(point,pointDimension,0,1);
			printf("%i------------- Punkt enthalten ------------:\n",(*ieq));	
			
			#endif
			return 1;
		}
	}
	#ifdef _DEBUG
	printf("------------- Ende der Ueberprfung von Punkt ------------:\n");
	printar(point,pointDimension,0,1);
	printf("------------- Punkt nicht enthalten ------------:\n");
	#endif
		
    return 0;
}  

SEXP SIneq(RAT * ieq, int noVars, int noIeq, int numerator){
	SEXP ans;
//	PROTECT(ans=allocVector(VECSXP,2));
//	PROTECT(ans=allocMatrix(REALSXP,noIeq,noVars+1));
	PROTECT(ans=allocMatrix(INTSXP,noIeq,noVars+1));
//	PROTECT(ansNum=allocMatrix(INTSXP,noIeq,noVars+1));
//	PROTECT(ansDen=allocMatrix(INTSXP,noIeq,noVars+1));	
	int i , j;	
	int column = noVars + 1;	
	for (i = 0; i < noIeq; i++){
		for (j = 0; j<column; j++) {    			
//			REAL(ans)[j*noIeq+i] = (ieq[(i)*column+j].num==0)?0:(double)ieq[(i)*column+j].num/ieq[(i)*column+j].den.i;			
			if (numerator) INTEGER(ans)[j*noIeq+i] = ieq[(i)*column+j].num; else INTEGER(ans)[j*noIeq+i] = ieq[(i)*column+j].den.i;						
		}				
	}
//	SET_VECTOR_ELT(ans,0,ansNum);
//	SET_VECTOR_ELT(ans,1,ansDen);
	UNPROTECT(1);
	return ans;
}


// ##############################################################################
// ##########################################################################TEST

void coneReduction(RAT * RSourceMatrix, int * dimRMatrix, RAT * RPointsToTest, int * RdimTestPoints, int * failurePoints, int ** matrixNumLeft, int ** matrixDenLeft, 
int ** matrixNumEliminated, int ** matrixDenEliminated, int ** oldRowsLeft, int ** oldRowsEliminated, int  ** dimLundE, char *** ieQch, int ** ieQ, SEXP * conesVar,SEXP * conesVarDen){
	/*	
	* @param  *pmatrix       ^= pointer to an double array (matrix) which stores the points (erfolgs/versagen) (row major) 
	*         *dimMatrix     ^= pointer to an int array which stores the dimension of pmatrix  (dimension of dimMatrix should be 2)
	*         *pointsToTest  ^= points to an double array which stores the points that will be tested 
	*         *dimPTT        ^= same as dimMatrix but for pointsToTest
	*         *failurePoints ^= failurePoints an array which stores for each point an 1 if it is an erfolgspunkt, 0 if it is a failurepoint
	* */
	
	
	
	int i,j;
	int dimMatrix  [2];
	RAT * matrix=0;
	
	dimMatrix[0] = dimRMatrix[0];
	dimMatrix[1] = dimRMatrix[1];
	
	matrix = (RAT *) RATallo(matrix,0,dimRMatrix[0]*dimRMatrix[1]);
		
	
    #ifdef _DEBUG//------------------------------------------------------------------------------------------
	printf("RMatrix:  \n");	
	printMatrix2v(RSourceMatrix,dimRMatrix);	
	printf("Matrix wird transponiert \n");
	#endif //------------------------------------------------------------------------------------------#
	
	for (i = 0; i <dimMatrix[0]; i++) {// wurde vorher tranponiert eingelesen
		for (j = 0; j <dimMatrix[1]; j++){
			matrix[i*dimMatrix[1]+j].num   = RSourceMatrix[i*dimMatrix[1]+j].num;
			matrix[i*dimMatrix[1]+j].den.i = RSourceMatrix[i*dimMatrix[1]+j].den.i;
			//I_RAT_Kuerzen(&(matrix[i*dimMatrix[1]+j]));
		}
	}	
	
	#ifdef _DEBUG    //------------------------------------------------------------------------------------------
	printar(matrix,dimMatrix[1],0,dimMatrix[0]);		
	printf("ausgefuehrt ---------\n");
	#endif //------------------------------------------------------------------------------------------#
	
	//versagenPunkte = versagPunkte(dimMatrix,failurePoints,matrix,&dimVersagen[0],&dimVersagen[1]);
	int dimFail[2];
	int dimSucc[2];
	RAT * matFail;
	RAT * matSucc;
	dimFail[1] = dimSucc[1] = dimMatrix[1];
	classifyPoints(matrix,dimMatrix,failurePoints,&matFail,&matSucc,&dimFail[0],&dimSucc[0]);
	#ifdef _DEBUG //------------------------------------------------------------------------------------------
	printf("es folgen die Versagenspunkte: \n");	
	printMatrix2v(matFail,dimFail);	
	printf("\n");	
	printf("es folgen die Erfolgspunkte: \n");	
	printMatrix2v(matSucc,dimSucc);	
	printf("\n");	
	printf("es folgen die PolyCone: \n");
    #endif //------------------------------------------------------------------------------------------#
	
	RAT * polyCones;
	polyCones=0;
	
	
	constructPolyCones(matSucc,dimSucc,matFail,dimFail,&polyCones);
	#ifdef _DEBUG //------------------------------------------------------------------------------------------
	for (i = 0; i < dimFail[0]; i++){
		printf("%i-ter PolyCone: \n",i);
		printar(&(polyCones[i * (dimSucc[0]*dimSucc[1])]),dimSucc[1],0,dimSucc[0]);
	}
	
	
	printf("\n");	
	printf("Die Test Punkte werden Transponiert und in RAT's konvertiert\n");
	#endif //------------------------------------------------------------------------------------------#
	
	
	RAT * pointsToTest = 0;
	pointsToTest = (RAT *) RATallo(RP pointsToTest,0,RdimTestPoints[0]*RdimTestPoints[1]);
	int dimPointsToTest[] = {RdimTestPoints[0],RdimTestPoints[1]};
	
	for (i = 0; i < dimPointsToTest[0]; i++){
		for (j = 0; j < dimPointsToTest[1]; j++){
			pointsToTest[i * dimPointsToTest[1] + j].num = 
				        RPointsToTest[i * dimPointsToTest[1] + j].num;
			pointsToTest[i * dimPointsToTest[1] + j].den.i = 
						RPointsToTest[i * dimPointsToTest[1] + j].den.i;
			//I_RAT_Kuerzen(&pointsToTest[i * dimPointsToTest[1] + j]);
		}
	}
	
	
	#ifdef _DEBUG//------------------------------------------------
	printar(pointsToTest,dimPointsToTest[1],0,dimPointsToTest[0]);	
	printf("\n--- ausgefuehrt \n");	
		
	printf("Konvertiere die Kegel in das Gleichungsformat: \n");
	#endif //----------------------------------------------------#
	RAT * ineqCones[dimFail[0]];
	int noOfFaces[dimFail[0]];
	for (i = 0; i < dimFail[0]; i++){
		
		int tempNoOfFaces;
		#ifdef _DEBUG //----------------------------------------
		printf("vor %iter Konvertierung \n",i);
	    #endif//-----------------------------------------------#
	    
    
	    
	    //''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
	    
		ineqCones[i] = 0;
		ineqCones[i] = testPorta(&polyCones[i * (dimSucc[0]*dimSucc[1])],dimSucc[0],dimSucc[1],&tempNoOfFaces);
		#ifdef _DEBUG  //----------------------------------
		printf("nach %iter Konvertierung \n",i);
	    #endif //------------------------------------------#		
		
		//''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
		
		
		
		
		(noOfFaces[i]) = tempNoOfFaces;
		#ifdef _DEBUG //------------------------------------------------------------------------------------------
		printf("%i-ter PolyCone: \n",i);
		printIneq((RAT*)ineqCones[i],dimFail[1],(int)noOfFaces[i]);
		printf("\n");//------------------------------------------------------------------------------------------#
		#endif
	}
	
	#ifdef _DEBUG //------------------------------------------------------------------------------------------
	printf("Die Gleichungen der Kegel lauten: \n");
	for (i = 0; i < dimFail[0]; i ++){
		printf("%i-ter PolyCone: \n",i);
		printIneq((RAT*)ineqCones[i],dimFail[1],(int)noOfFaces[i]);
		printf("\n");
	}
	printf("\n");
	
	printf("Wandle Failure Points ins RAT format\n");
	#endif//------------------------------------------------------------------------------------------#
	//RAT * ratFailurePoints;
	
	//ratFailurePoints = (RAT*) RATallo(ratFailurePoints,0,dimFail[0]*dimFail[1]);
	#ifdef _DEBUG 
	printf(" \n\n*******************#######################-----------------------***********************");
	#endif		
	//for (i = 0; i < dimFail[0]; i++){
	//	for (j = 0; j < dimFail[1]; j++){
	//		ratFailurePoints[i * dimFail[1] + j].num = 
	//		loRound(matFail[i* dimFail[1] + j]* precision);
	//		ratFailurePoints[i * dimFail[1] + j].den.i = 
	//			precision;
	//		I_RAT_Kuerzen(&ratFailurePoints[i * dimFail[1] + j]);
	//	}
	//}
	 
	#ifdef _DEBUG_TODO 
	printf(" \n\n**********************************************\nDie enthaltene Punkte sind:   \n");
	for (i = 0; i < dimPointsToTest[0]; i++){
		if (testPointInUnionOfPolyCones(ineqCones,matFail,noOfFaces,				// Bei konvertierung in CPP Datei hab ich hier die Sternchen 
										dimFail[0],&pointsToTest[i*dimPointsToTest[1]],		// bei dem 1. und 3. Argument ergaenzt
										dimPointsToTest[1],&ieq)){
			printar(&pointsToTest[i*dimPointsToTest[1]],dimPointsToTest[1],0,1);
			printf(" enthalten \n");
		} else {
			printar(&pointsToTest[i*dimPointsToTest[1]],dimPointsToTest[1],0,1);
			printf(" nicht enthalten \n");
		}		
	} 
	#endif		
	//'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''#
	int left = 0;
	int	ieq;	//''''''''''''''''''''''''' 
	int mem_TestPointInUnionOfPolyCones_Result[dimPointsToTest[0]];
	for (i = 0; i < dimPointsToTest[0]; i++){
		mem_TestPointInUnionOfPolyCones_Result[i] = 
			testPointInUnionOfPolyCones(ineqCones,matFail,noOfFaces,				// Bei konvertierung in CPP Datei hab ich hier die Sternchen 
										dimFail[0],&pointsToTest[i*dimPointsToTest[1]],		// bei dem 1. und 3. Argument ergaenzt
										dimPointsToTest[1] , &ieq) == 0;
		if (mem_TestPointInUnionOfPolyCones_Result[i]){
				left++;			 
		}
	}
	
		
	(*dimLundE)       = (int *) malloc (3 * sizeof(int));	    
    (*dimLundE)[0]   = left;
	(*dimLundE)[1]   = dimPointsToTest[0] - left ;
	(*dimLundE)[2]   = dimFail[0];

    
    (*oldRowsLeft)   = (int *) malloc (left * sizeof(int));
	(*matrixNumLeft) = (int *) malloc (left * dimPointsToTest[1] * sizeof(int));
	(*matrixDenLeft) = (int *) malloc (left * dimPointsToTest[1] * sizeof(int));
	       	
	RAT * par2;		
	int  temp2;
	int  temp3=0;
	
	
	for (i = 0; i < dimPointsToTest[0]; i++){
		if (mem_TestPointInUnionOfPolyCones_Result[i]){
		    (*oldRowsLeft)[temp3] = i;     
        	par2 = &pointsToTest[i*dimPointsToTest[1]];    
			for (temp2 = 0; temp2 < dimPointsToTest[1]; temp2++){			   	
				(*matrixNumLeft)[temp3*dimPointsToTest[1] + temp2] = (int) par2[temp2].num;				
				(*matrixDenLeft)[temp3*dimPointsToTest[1] + temp2] = (int) par2[temp2].den.i;								
				#ifdef DEBUG_R_RUECKGABE
				printf("%i/%i ", (int) par2[temp2].num,(int) par2[temp2].den.i);
				#endif
			}
			#ifdef DEBUG_R_RUECKGABE
			printf("\n");
			#endif
			temp3++;	
		}	
		
	}
	#ifdef DEBUG_R_RUECKGABE
	printf("\nleft: %i, temp3: %i\n", left,temp3);
	#endif
	temp3 = 0;
	
    int eliminated = dimPointsToTest[0] - left;
    (*oldRowsEliminated)  = (int *) malloc (eliminated * sizeof(int));
    (*matrixNumEliminated)= (int *) malloc (eliminated * dimPointsToTest[1] * sizeof(int));
	(*matrixDenEliminated)= (int *) malloc (eliminated * dimPointsToTest[1] * sizeof(int));
	#ifdef DEBUG_R_RUECKGABE
	printf("Eliminated");
	#endif
	for (i = 0; i < dimPointsToTest[0]; i++){
		if (mem_TestPointInUnionOfPolyCones_Result[i]==0){
		    (*oldRowsEliminated)[temp3] = i;     
        	par2 = &pointsToTest[i*dimPointsToTest[1]];
    
			for (temp2 = 0; temp2 < dimPointsToTest[1]; temp2++){			   	
				(*matrixNumEliminated)[temp3*dimPointsToTest[1] + temp2] = (int) par2[temp2].num;				
				(*matrixDenEliminated)[temp3*dimPointsToTest[1] + temp2] = (int) par2[temp2].den.i;
				#ifdef DEBUG_R_RUECKGABE
				printf("%i/%i ", (int) par2[temp2].num,(int) par2[temp2].den.i);
				#endif	
			}
			#ifdef DEBUG_R_RUECKGABE
			printf("\n");
			#endif
			temp3++;	
		}	
	}
	
	
	
	// ieQ[i] enthaelt die kegelNummer, die grund fuer die Elimination der i.te Knoten war   
    (*ieQ) = (int *) malloc (eliminated * sizeof(int));
    int temp4 = 0;
    for (i = 0; i < dimPointsToTest[0]; ++i) {
    	ieq = 0;
		if (testPointInUnionOfPolyCones(ineqCones,matFail,noOfFaces, dimFail[0],&pointsToTest[i*dimPointsToTest[1]], dimPointsToTest[1], &ieq) ){			
			 (*ieQ)[temp4++] = ieq;	
		}										
	}
			
	#ifdef _DEBUG 
	for (i = 0; i < eliminated; ++i) {
		printf("\nindex---------------->>>>>>>%i  ", (*ieQ)[i]);		
	}
	#endif
	
	
	PROTECT((*conesVar) = allocVector(VECSXP,dimFail[0]));
	PROTECT((*conesVarDen) = allocVector(VECSXP,dimFail[0]));

	
 	(*ieQch)  = (char **) malloc (eliminated * sizeof(char*));

	for (i = 0; i < dimFail[0]; ++i) {		
		
		SET_VECTOR_ELT((*conesVar),i,SIneq((RAT*)ineqCones[i],dimFail[1],(int)noOfFaces[i],1));
		SET_VECTOR_ELT((*conesVarDen),i,SIneq((RAT*)ineqCones[i],dimFail[1],(int)noOfFaces[i],0));		
	 	(*ieQch)[i] = stringIneq((RAT*)ineqCones[i],dimFail[1],(int)noOfFaces[i]);
			
	}

	#ifdef _DEBUG 
	for (i = 0; i < dimFail[0]; ++i) {
		printf((*ieQch)[i]);		 			
	}
	#endif 
		
}








/*
	
#ifndef COMPILE_AS_LIB
int main( int argc, char *argv[] ){
	int m = 1;	
	int RdimTestPoints[] = {3,216};
	double RPointsToTest[RdimTestPoints[1]* RdimTestPoints[0]]; 
	
	double  fac1[] = { 0.8 , 0.4 , 0.2 , 0.2 , 1.0 , 0.2 , 0.8, 0.2 };//, 0,0};		
	double  fac2[] = { 0.2 , 0.2 , 0.2 , 0.8 , 0.2 , 0.2 ,  0.8,  0.8 };//, 0 ,0};
	double  fac3[] = { 0.8 , 0.0 , 1.0 , 0.8 , 0.2 , 0.4 , 0.2,  0.2  };//, 0 ,0};	
	//int failurePoints[] = { 0,1,0,1,0,0,1,0,1,0}; //, 0, 0};//, 0,   1};
	//int failurePoints[] = { 0,0,0,1,1,0,0,0};
	int failurePoints[] = { 0,1,0,0,1,0,0,0};
	int dimRMatrix  [2] =  {3,8};//{ 3, 5};
	int i;
	
	double x,y,z;
	int pC = 0;
	for (x = 0; x <= 1.0; x+=0.2){
		for (y = 0; y <= 1.0; y+= 0.2){
			for (z = 0; z <= 1.0; z+= 0.2){
				RPointsToTest[RdimTestPoints[1] * 0 + pC] = x;
				RPointsToTest[RdimTestPoints[1] * 1 + pC] = y;
				RPointsToTest[RdimTestPoints[1] * 2 + pC] = z;
				pC++;
			}
		}
	}
	
	double RSourceMatrix [dimRMatrix[0]*dimRMatrix[1]];
	for (i = 0; i < dimRMatrix[1]; i++){
		RSourceMatrix[i+dimRMatrix[1]*0] = fac1[i];
	}
	for (i = 0; i < dimRMatrix[1]; i++){
		RSourceMatrix[i+dimRMatrix[1]*1] = fac2[i];
	}
	for (i = 0; i < dimRMatrix[1]; i++){
		RSourceMatrix[i+dimRMatrix[1]*2] = fac3[i];
	}	
	
	int * dimANS = 0;	
	int * matrixNum = 0;
	int * matrixDen = 0;
	
	coneReduction((double*)&RSourceMatrix,(int*)&dimRMatrix,(double*)&RPointsToTest,(int*)&RdimTestPoints,failurePoints,
					&dimANS,&matrixNum,&matrixDen);
	int j;
	
	for (i = 0; i < dimANS[0]; i++){
		for (j = 0; j < dimANS[1]; j++){
			printf("%i/%i ",matrixNum[i*dimANS[1]+j],matrixDen[i*dimANS[1]+j]);
		}
		printf("\n");
	}
	
	return 0;
}
#endif
*/
//#ifdef COMPILE_AS_LIB 
 
int readRAT(const char * pRatStr, RAT * pInRat){
	char * ptr;
	char ch;
	(*pInRat).num = (int) strtol(pRatStr,&ptr,10);
	if (ptr == pRatStr){
		#ifdef _DEBUG   
		printf("invalid format");
		#endif 
		return 0;
	}
	/* position in_line after the number */
	pRatStr = ptr;
	do 
	{
		ch = *(pRatStr++);
	} 
	while (ch == ' ' || ch == '\t');
	if (ch  == '/') 
	{
		(*pInRat).den.i = (int) strtol(pRatStr,&ptr,10);
	    if (ptr == pRatStr || (*pInRat).den.i <= 0){
	    	#ifdef _DEBUG
			printf("invalid format or negative denumerator");
			#endif 
			return 0;
		}
	} else {
		/* default denominator = 1 */
	    (*pInRat).den.i = 1;
	}
	return 1;
}


SEXP rConeReduction(SEXP pRSourceMatrix, SEXP pRPointsToTest, SEXP pfailurePoints){
	int i,j;
	// Convert SourceMatrix
	int dimRMatrix[2];
		
	// Wenn die Matrix Transponiert eingelesen werden soll dann so :
	//dimRMatrix[0] = INTEGER(getAttrib(pRSourceMatrix,R_DimSymbol))[0];
	//dimRMatrix[1] = INTEGER(getAttrib(pRSourceMatrix,R_DimSymbol))[1];	
			
	//for(i=0;i<dimRMatrix[0];i++) {							
	//	for(j=0;j<dimRMatrix[1];j++) {
	//		if (isReal(pRSourceMatrix)==1) {
	//		   RSourceMatrix[i*dimRMatrix[1]+j].num=
	//			   loRound(REAL(pRSourceMatrix)[i+j*dimRMatrix[0]] * precision);
	//		   RSourceMatrix[i*dimRMatrix[1]+j].den.i = precision;
	//		} else if (isInteger(pSourceMatrix) == 1) {
	//		   RSourceMatrix[i*dimRMatrix[1]+j].num=
	//			   loRound(((double)INTEGER(pRSourceMatrix)[i+j*dimRMatrix[0]]) *precision);	
	//		   RSourceMatrix[i*dimRMatrix[1]+j].den.i = precision;
	//	   	}
	//		I_RAT_Kuerzen(&RSourceMatrix[i*dimRMatrix[1]+j]);
	//	}
	//}  

	 
	dimRMatrix[0] = INTEGER(getAttrib(pRSourceMatrix,R_DimSymbol))[0]; // evt. noch umdrehen
	dimRMatrix[1] = INTEGER(getAttrib(pRSourceMatrix,R_DimSymbol))[1];	
	
	RAT RSourceMatrix [dimRMatrix[0]*dimRMatrix[1]];		
	
	#ifdef _DEBUG
	printf("dimRMat: %i,%i",dimRMatrix[0],dimRMatrix[1]);
	#endif
	for(i=0;i<dimRMatrix[0];i++) {							
		for(j=0;j<dimRMatrix[1];j++) {
			#ifdef _DEBUG
			printf ("[%i,%i]",i,j);
			#endif
			if (isReal(pRSourceMatrix)==1) {
			   RSourceMatrix[i*dimRMatrix[1]+j].num= // wichtig R speichert matrizen spaltenweise 
				   loRound((REAL(pRSourceMatrix)[i + j*dimRMatrix[0]]) * (double)precision);
			   RSourceMatrix[i*dimRMatrix[1]+j].den.i = precision;
			} else if (isInteger(pRSourceMatrix)) {
			   RSourceMatrix[i*dimRMatrix[1]+j].num=
				   loRound(((double)INTEGER(pRSourceMatrix)[i+j*dimRMatrix[0]]) * (double)precision);	
			   RSourceMatrix[i*dimRMatrix[1]+j].den.i = precision;
		   	} else if (isString(pRSourceMatrix)){
		   		const char * copy;
		   		PROTECT(pRSourceMatrix);// = AS_CHARACTER(pRSourceMatrix));
		   		// allocate memory to Pmychar[0], Pmychar[1]:
		   		//copy = R_alloc(strlen(CHAR(STRING_ELT(pRSourceMatrix, i+j*dimRMatrix[0]))), sizeof(char));
		   		copy = CHAR(STRING_ELT(pRSourceMatrix, i+j*dimRMatrix[0]));
		   		
		   		#ifdef _DEBUG
		   		printf("copy: %s",copy);
		   		#endif 
		   		//strcpy(copy, CHAR(STRING_ELT(pRSourceMatrix, i+j*dimRMatrix[0])));
		   		readRAT(copy ,&RSourceMatrix[i*dimRMatrix[1]+j]);
		   		
		   		UNPROTECT(1);

		   	}
			I_RAT_Kuerzen(&RSourceMatrix[i*dimRMatrix[1]+j]);
		}
	} 
	
	
	// Convert PointsToTest	
	int RdimTestPoints[2];	
	
//	RdimTestPoints[0] = INTEGER(getAttrib(pRPointsToTest,R_DimSymbol))[0];
//	RdimTestPoints[1] = INTEGER(getAttrib(pRPointsToTest,R_DimSymbol))[1];	
	
//	RAT RPointsToTest [RdimTestPoints[0]*RdimTestPoints[1]];		
//	for(i=0;i<RdimTestPoints[0];i++) {							
//		for(j=0;j<RdimTestPoints[1];j++) {
//			if (isReal(pRPointsToTest)==1){ 
//				RPointsToTest[i*RdimTestPoints[1]+j].num=
//					loRound(REAL(pRPointsToTest)[i+RdimTestPoints[0]*j] *precision);
//				RPointsToTest[i*RdimTestPoints[1]+j].den.i = precision;			
//			}else{
//				RPointsToTest[i*RdimTestPoints[1]+j].num=
//					loRound(((double)INTEGER(pRPointsToTest)[i+RdimTestPoints[0]*j])*precision);
//				RPointsToTest[i*RdimTestPoints[1]+j].den.i = precision;								
//			}
//			I_RAT_Kuerzen(&RPointsToTest[i*RdimTestPoints[1]+j]);
//		}
//	}	
	
	RdimTestPoints[0] = INTEGER(getAttrib(pRPointsToTest,R_DimSymbol))[0]; // evtl noch umdrehen
	RdimTestPoints[1] = INTEGER(getAttrib(pRPointsToTest,R_DimSymbol))[1];	
	
	RAT RPointsToTest [RdimTestPoints[0]*RdimTestPoints[1]];		
	for(i=0;i<RdimTestPoints[0];i++) {							
		for(j=0;j<RdimTestPoints[1];j++) {
			if (isReal(pRPointsToTest)==1){ 
				RPointsToTest[i*RdimTestPoints[1]+j].num=
					loRound(REAL(pRPointsToTest)[j*RdimTestPoints[0]+i] * (double)precision); // wichtig R speichert matrix spaltenweise
				RPointsToTest[i*RdimTestPoints[1]+j].den.i = precision;			
			}else if (isInteger(pRSourceMatrix)) {
				RPointsToTest[i*RdimTestPoints[1]+j].num=
					loRound(((double)INTEGER(pRPointsToTest)[j*RdimTestPoints[0]+i])* (double)precision);
				RPointsToTest[i*RdimTestPoints[1]+j].den.i = precision;								
			} else if (isString(pRSourceMatrix)){
		   		const char * copy;
		   		PROTECT(pRPointsToTest);// = AS_CHARACTER(pRSourceMatrix));
		   		// allocate memory to Pmychar[0], Pmychar[1]:
		   		//copy = R_alloc(strlen(CHAR(STRING_ELT(pRSourceMatrix, i+j*dimRMatrix[0]))), sizeof(char));
		   		copy = CHAR(STRING_ELT(pRPointsToTest, j*RdimTestPoints[0]+i));
		   		#ifdef _DEBUG
		   		printf("copy: %s",copy);
		   		#endif 
		   		//strcpy(copy, CHAR(STRING_ELT(pRSourceMatrix, i+j*dimRMatrix[0])));
		   		readRAT(copy ,&RPointsToTest[i*RdimTestPoints[1]+j]);
		   		
		   		UNPROTECT(1);

		   	}
			I_RAT_Kuerzen(&RPointsToTest[i*RdimTestPoints[1]+j]);
		}
	}  
	

		
	// Convert failurePoints	
	int failDim=length(pfailurePoints);	
	int failurePoints[failDim];	
	for (i = 0; i < failDim; i++){
		failurePoints[i] = INTEGER(pfailurePoints)[i];
	}	
	
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	
	SEXP ans,matrixLEFT,dimNames,dimNames2,rowNames,rowNames2,colNames,colNames2,matrixELIMINATED,cones; 
		
	SEXP conesVar,conesVarDen;	
		
	char  ** ieQch;
	
	int * ieQ; 
	
	int * dimLundE;	
		
	int * matrixNumLeft = 0;
	int * matrixDenLeft = 0;
	int * oldRowsLeft   = 0;

    int * matrixNumEliminated = 0; 
	int * matrixDenEliminated = 0;
	int * oldRowsEliminated   = 0;	
		
	coneReduction(RSourceMatrix,dimRMatrix,RPointsToTest,RdimTestPoints,failurePoints, &matrixNumLeft, &matrixDenLeft, 
	&matrixNumEliminated,&matrixDenEliminated, &oldRowsLeft, &oldRowsEliminated, &dimLundE, &ieQch, &ieQ,&conesVar,&conesVarDen);
	
	int dim           	= RdimTestPoints[1];
	int dimLeft       	= dimLundE[0];
	int dimEliminated 	= dimLundE[1];
	int dimFail			= dimLundE[2];
	
	
	PROTECT(ans=allocVector(VECSXP,4));
	
		// left candidates
    PROTECT(matrixLEFT = allocMatrix(REALSXP,dimLeft,dim));
    PROTECT(dimNames = allocVector(VECSXP,2));
	PROTECT(rowNames = allocVector(INTSXP,dimLeft));
	PROTECT(colNames = allocVector(INTSXP,dim));	    
	int t1,t2;
//	 	dimNames=getAttrib(pRSourceMatrix,R_DimNamesSymbol);
	for(t1=0;t1<dimLeft;t1++) INTEGER(rowNames)[t1]=oldRowsLeft[t1]+1;							
	for(t2=0;t2<dim;t2++) INTEGER(colNames)[t2]=t2+1;
	SET_VECTOR_ELT(dimNames,0,rowNames);
	SET_VECTOR_ELT(dimNames,1,colNames);
	dimnamesgets(matrixLEFT,dimNames);

		#ifdef DEBUG_R_RUECKGABE
	 	printf("RConeReduction: matrixNumLeft / matrixDenLeft: \n")
	 	for(t1=0;t1<dimLeft;t1++) {							
			for(t2=0;t2<dim;t2++) {	 	
				printf("%i/%i ",(int)matrixNumLeft[t1*dim+t2],(int)matrixDenLeft[t1*dim+t2]);
			}
			printf("\n");
	 	}		  
		#endif
	 	
	 	for(t1=0;t1<dimLeft;t1++) {							
			for(t2=0;t2<dim;t2++) {	 	
			REAL(matrixLEFT)[t2*dimLeft+t1] = (double)matrixNumLeft[t1*dim+t2]/matrixDenLeft[t1*dim+t2] ;
			}
	 	}		  
			
		// excluded candidates	
	    PROTECT(matrixELIMINATED = allocMatrix(REALSXP,dimEliminated,dim+1));
	    PROTECT(dimNames2 = allocVector(VECSXP,2));	    
	    PROTECT(rowNames2 = allocVector(INTSXP,dimEliminated));
	    PROTECT(colNames2 = allocVector(INTSXP,dim+1));	    	    
	 	for(t1=0;t1<dimEliminated;t1++) INTEGER(rowNames2)[t1]=oldRowsEliminated[t1]+1;							
	 	for(t2=0;t2<dim+1;t2++) INTEGER(colNames2)[t2]=t2+1;
	 	SET_VECTOR_ELT(dimNames2,0,rowNames2);
	 	SET_VECTOR_ELT(dimNames2,1,colNames2);	 	
	 	dimnamesgets(matrixELIMINATED,dimNames2);
	 	for(t1=0;t1<dimEliminated;t1++) {							
			for(t2=0;t2<dim;t2++) {	 	
			REAL(matrixELIMINATED)[t2*dimEliminated+t1] = (double)matrixNumEliminated[t1*dim+t2]/matrixDenEliminated[t1*dim+t2] ;
			//REAL(matrixELIMINATED)[t1*dim+t2] = (double)matrixNumEliminated[t1*dim+t2]/matrixDenEliminated[t1*dim+t2] ;
			}
	 	}		  
	 	for(t1=0;t1<dimEliminated;t1++) {							
			REAL(matrixELIMINATED)[dim*dimEliminated+t1] = ieQ[t1]+1 ;
			
			INTEGER(rowNames2)[t1]=ieQ[t1]+1;
			
	 	}		  

		// cones
		
		PROTECT(cones=allocVector(STRSXP,dimFail));
	 	for(t1=0;t1<dimFail;t1++) SET_STRING_ELT(cones,t1,mkChar(ieQch[t1]));
		
		 	
	 	SET_VECTOR_ELT(ans,0,rowNames);	 		
	 	SET_VECTOR_ELT(ans,1,rowNames2);	 		
	 	SET_VECTOR_ELT(ans,2,conesVar);	 			 			 		
	 	SET_VECTOR_ELT(ans,3,conesVarDen);	 	
	 		 
	UNPROTECT(12);
	
	return ans;		
	
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
}
 

SEXP getMaxINT(){	
	SEXP maxINT;
	PROTECT(maxINT = allocVector(INTSXP,1));
	int max = INT_MAX;
	INTEGER (maxINT)[0] = max;	
	UNPROTECT(1);
	return maxINT;
}
 
 
	    	     
SEXP rConeReductionINT(SEXP pRSourceMatrixNum,SEXP pRSourceMatrixDen, SEXP pRPointsToTestNum, SEXP pRPointsToTestDen, SEXP pfailurePoints){
	int i,j;
	// Convert SourceMatrix ,which consist of two integer Matrix 
	int dimRSMatrix[2];
	dimRSMatrix[0] = INTEGER(getAttrib(pRSourceMatrixDen,R_DimSymbol))[0];
	dimRSMatrix[1] = INTEGER(getAttrib(pRSourceMatrixDen,R_DimSymbol))[1];
			
	RAT * RSourceMatrix = 0;
	RSourceMatrix = (RAT *) RATallo(RSourceMatrix,0,dimRSMatrix[0]*dimRSMatrix[1]);	
	for (i = 0; i < dimRSMatrix[0]; i++){
		for (j = 0; j < dimRSMatrix[1]; j++){
			RSourceMatrix[i * dimRSMatrix[1] + j].num   = (int) REAL(pRSourceMatrixNum)[i+ j*dimRSMatrix[0]];
			RSourceMatrix[i * dimRSMatrix[1] + j].den.i = (int) REAL(pRSourceMatrixDen)[i+ j*dimRSMatrix[0]];		
			I_RAT_Kuerzen(&RSourceMatrix[i * dimRSMatrix[1] + j]);
		}
	}

	// Convert PointsToTest		
	int RdimTestPoints[2];	
	RdimTestPoints[0] = INTEGER(getAttrib(pRPointsToTestDen,R_DimSymbol))[0];
	RdimTestPoints[1] = INTEGER(getAttrib(pRPointsToTestDen,R_DimSymbol))[1];
		
	
	RAT * pointsToTest = 0;
	pointsToTest = (RAT *) RATallo(pointsToTest,0,RdimTestPoints[0]*RdimTestPoints[1]);
	//int dimPointsToTest[] = {RdimTestPoints[1],RdimTestPoints[0]};

    for (i = 0; i < RdimTestPoints[0]; i++){
		for (j = 0; j < RdimTestPoints[1]; j++){
			pointsToTest[i * RdimTestPoints[1] + j].num     = (int) REAL(pRPointsToTestNum)[i + j*RdimTestPoints[0]];
			pointsToTest[i * RdimTestPoints[1] + j].den.i   = (int) REAL(pRPointsToTestDen)[i + j*RdimTestPoints[0]];
			I_RAT_Kuerzen(&pointsToTest[i * dimRSMatrix[1] + j]);
		}
	}

	 	
	// Convert failurePoints	
	int failDim=length(pfailurePoints);	
	int failurePoints[failDim];	
	for (i = 0; i < failDim; i++){
		failurePoints[i] = INTEGER(pfailurePoints)[i];
	}	
	
	
	
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	
	
	SEXP ans,matrixLeftNum,matrixLeftDen,dimNames,dimNames2,rowNames,rowNames2,colNames,colNames2,matrixEliminatedNum,matrixEliminatedDen,cones; 
		
	SEXP conesVar,conesVarDen;	
		
	char  ** ieQch;
	
	int * ieQ; 
	
	int * dimLundE;	
		
	int * matrixNumLeft = 0;
	int * matrixDenLeft = 0;
	int * oldRowsLeft   = 0;

    int * matrixNumEliminated = 0; 
	int * matrixDenEliminated = 0;
	int * oldRowsEliminated   = 0;	

		
	coneReduction(RSourceMatrix,dimRSMatrix,pointsToTest,RdimTestPoints,failurePoints, &matrixNumLeft, &matrixDenLeft, 
	&matrixNumEliminated,&matrixDenEliminated, &oldRowsLeft, &oldRowsEliminated, &dimLundE, &ieQch, &ieQ,&conesVar,&conesVarDen);
	
	
	int dim           	= RdimTestPoints[1];
	int dimLeft       	= dimLundE[0];
	int dimEliminated 	= dimLundE[1];
	int dimFail			= dimLundE[2];
	
	
	PROTECT(ans=allocVector(VECSXP,4));
	
		// left candidates
	    PROTECT(matrixLeftNum = allocMatrix(INTSXP,dimLeft,dim));
	    PROTECT(matrixLeftDen = allocMatrix(INTSXP,dimLeft,dim));
	    PROTECT(dimNames = allocVector(VECSXP,2));
	    PROTECT(rowNames = allocVector(INTSXP,dimLeft));
	    PROTECT(colNames = allocVector(INTSXP,dim));	    
	 	int t1,t2;
//	 	dimNames=getAttrib(pRSourceMatrix,R_DimNamesSymbol);
	 	for(t1=0;t1<dimLeft;t1++) INTEGER(rowNames)[t1]=oldRowsLeft[t1]+1;							
	 	for(t2=0;t2<dim;t2++) INTEGER(colNames)[t2]=t2+1;
	 	SET_VECTOR_ELT(dimNames,0,rowNames);
	 	SET_VECTOR_ELT(dimNames,1,colNames);
	 	dimnamesgets(matrixLeftNum,dimNames);
	 	dimnamesgets(matrixLeftDen,dimNames);
	 	for(t1=0;t1<dimLeft;t1++) {							
			for(t2=0;t2<dim;t2++) {	 	
			INTEGER(matrixLeftNum)[t2*dimLeft+t1] = matrixNumLeft[t1*dim+t2];
			INTEGER(matrixLeftDen)[t2*dimLeft+t1] = matrixDenLeft[t1*dim+t2];
			}
	 	}		  
	
			
		// excluded candidates	 
	    PROTECT(matrixEliminatedNum = allocMatrix(INTSXP,dimEliminated,dim+1));
	    PROTECT(matrixEliminatedDen = allocMatrix(INTSXP,dimEliminated,dim+1));
	    PROTECT(dimNames2 = allocVector(VECSXP,2));	    
	    PROTECT(rowNames2 = allocVector(INTSXP,dimEliminated));
	    PROTECT(colNames2 = allocVector(INTSXP,dim+1));	    	    
	 	for(t1=0;t1<dimEliminated;t1++) INTEGER(rowNames2)[t1]=oldRowsEliminated[t1]+1;							
	 	for(t2=0;t2<dim+1;t2++) INTEGER(colNames2)[t2]=t2+1;
	 	SET_VECTOR_ELT(dimNames2,0,rowNames2);
	 	SET_VECTOR_ELT(dimNames2,1,colNames2);	 	
	 	dimnamesgets(matrixEliminatedNum,dimNames2);
	 	dimnamesgets(matrixEliminatedDen,dimNames2);
	 	for(t1=0;t1<dimEliminated;t1++) {							
			for(t2=0;t2<dim;t2++) {	 	
			INTEGER(matrixEliminatedNum)[t2*dimEliminated+t1] = matrixNumEliminated[t1*dim+t2];
			INTEGER(matrixEliminatedDen)[t2*dimEliminated+t1] = matrixDenEliminated[t1*dim+t2];
			}
	 	}		  
	 	for(t1=0;t1<dimEliminated;t1++) {							
			INTEGER(matrixEliminatedNum)[dim*dimEliminated+t1] = ieQ[t1]+1 ;
			INTEGER(matrixEliminatedDen)[dim*dimEliminated+t1] = ieQ[t1]+1 ;
			
			INTEGER(rowNames2)[t1]=ieQ[t1]+1;
	 	}		  

		// cones
		
		PROTECT(cones=allocVector(STRSXP,dimFail));
	 	for(t1=0;t1<dimFail;t1++) SET_STRING_ELT(cones,t1,mkChar(ieQch[t1]));
		
		 	
	 	SET_VECTOR_ELT(ans,0,rowNames);
	 	SET_VECTOR_ELT(ans,1,rowNames2);
	 	SET_VECTOR_ELT(ans,2,conesVar);	 			 			 		
	 	SET_VECTOR_ELT(ans,3,conesVarDen);	 	
	 		 
	UNPROTECT(14);
	
	return ans;		

	
	//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	
	
}
// TODO FileName2! isNull überprüfen etc pp
SEXP rPortaInterface(SEXP fileName, SEXP options,SEXP fileName2){
	
/*	printf("rPortaInterface betreten");
	int length = LENGTH(options);
		
	PROTECT(options = AS_CHARACTER(options));
	
	printf("laenge %i", length);
		
	
	char * temp_Options[length];
	int i = 0; 
	
	for (i = 0; i < length; i++){
		temp_Options[i] = R_alloc(strlen(CHAR(STRING_ELT(options, i))), sizeof(char));
		strcpy(temp_Options[i], CHAR(STRING_ELT(options, i)));
		printf("%s option %i \n",temp_Options[i],i);
	}
	
	char * porta_FileName;
	
	PROTECT(fileName = AS_CHARACTER(fileName));
	porta_FileName = R_alloc(strlen(CHAR(STRING_ELT(fileName, 0))), sizeof(char));
	strcpy(porta_FileName, CHAR(STRING_ELT(fileName, 0)));
	 
	printf("\nlaenge:%i\n",length);*/
	 
	int numberOfOptions=isNull(fileName2)?3:4;
	char * porta_Options[numberOfOptions];
	
	char * nothing = "-nothing";
	
	porta_Options[2] = (char *)CHAR(STRING_ELT(fileName,0));
	if (numberOfOptions==4) porta_Options[3] = (char *)CHAR(STRING_ELT(fileName2,0));
/*	porta_Options[2] = porta_FileName;
	char convert_option [length+1];
	convert_option[0] = '-';
	
	for (i = 0; i < length; i++){
		convert_option[i+1] = (char) (temp_Options[i])[1];
	}*/
	
	porta_Options[0] = nothing;
//	porta_Options[1] = convert_option;
	porta_Options[1] = (char *)CHAR(STRING_ELT(options,0));
	
	#ifdef _DEBUG 
	if (numberOfOptions==4) printf("\n FileName2: \"%s\"\n",porta_Options[3]);
	printf("\n FileName: \"%s\"\n",porta_Options[2]); 
	printf("\n Options: \"%s\"\n",porta_Options[1]);
	#endif 		
	
	if (numberOfOptions==4) {
		// Nicht die eleganteste aber eine sehr einfache Lösung:
		if (strcmp(porta_Options[2],porta_Options[3])==0) {
			//vint
			porta_Options[3]="";
			numberOfOptions=3;
		}
		validmain(numberOfOptions,porta_Options);
	} else {	
		portamain(numberOfOptions,porta_Options);
	}	
	
	UNPROTECT(1);
	 
	
	//printf("Debug: Rueckgabe von IEQPOI in POI.c \n");

	
	return IEQPOI;
}
 
 //----------------------- C Anbindung 
 /* 
R_CMethodDef cMethods[] = {
	{"coneReduction",&coneReduction,5,{REALSXP,INTSXP,REALSXP,INTSXP,LGLSXP}},
	{NULL,NULL,0}
}; // muss evtl. ausserhalb der Methode deklariert werden damit danach auch noch verfgbar <-- testen !  
 
 
 
 R_CallMethodDef callMethods[] = {
 	{"rConeReduction",&rConeReduction,5},
 	{NULL,NULL,0}
 };
 
 
 
void R_init_RPorta(DllInfo *info){
	
	// wird beim laden der DLL als erstes aufgerufen 
	// hier koennen die Funktionen registriert werden die nachher ueber 
	// r mittels .C / .Call ansprechbar sein sollen 
	 
	
	
	R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}
 
void R_unload_RPorta(DllInfo *info)
{
	
	// wird beim entladen der DLL aufgerufen 
	// hier koennen die Resourcen die waerend des Rechen vorgangs benutzt werden wieder
	// freigegeben werden 
	 
} 

void R_init_sfb475a5(DllInfo *info){
	
	// wird beim laden der DLL als erstes aufgerufen 
	// hier koennen die Funktionen registriert werden die nachher ueber 
	// r mittels .C / .Call ansprechbar sein sollen 
	 
	
	
	R_registerRoutines(info,NULL,callMethods,NULL,NULL);
}
 
void R_unload_sfb475a5(DllInfo *info)
{
	
	// wird beim entladen der DLL aufgerufen 
	// hier koennen die Resourcen die waerend des Rechen vorgangs benutzt werden wieder
	// freigegeben werden 
	 
} 
  */



//#endif
