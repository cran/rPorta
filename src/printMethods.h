/*
 * This file is part of RPorta. For licensing and copyright information
 * please see the file COPYING in the root directory of this
 * distribution or contact <robin.nunkesser@tu-dortmund.de>.
 */

#include "porta.h"
#include "common.h"
#include "string.h"
 

//#define PRINT_DEBUG

void printar(RAT * par, int column, int first, int last){
	int ie,i;
	for (ie = first; ie < last; ie++) 
    {
		for (i = 0; i < column; i++){
			printf("%i/%i \t" , (int) par[(ie)*column+i].num, par[(ie)*column+i].den.i);
			
		}
		printf("\n");
	}
}


void printIneq(RAT * ieq, int noVars, int noIeq){
	int i , j;
	int column = noVars + 1;
	// ieq a1 * x1 + a2 * x2 + ... + an * xn <= b
	for (i = 0; i < noIeq; i++){
		for (j = 0; j<noVars-1; j++) {
			printf("%i/%i*x%i + " , (int)ieq[(i)*column+j].num, (int)ieq[(i)*column+j].den.i,j);
		}
		printf("%i/%i*x%i <= " , (int)ieq[(i)*column+noVars-1].num, (int)ieq[(i)*column+noVars-1].den.i,noVars-1);
		printf("%i/%i \n" , (int)ieq[(i)*column+noVars].num, (int)ieq[(i)*column+noVars].den.i);
	}
}

	
 char * konkatenation(char * eins, char * zwei){                      
	    	int length = strlen(eins) + strlen(zwei);
			char * ktn =  malloc ((length+1 )* sizeof(char)); 
			
	    	int i ;
	    	for (i =0 ; i < strlen(eins); i++) 
				ktn[i] = eins[i]; 
			int temp =0;
	    	
			for (i = strlen(eins); i < length; i++) 
				ktn[i] = zwei[ temp++];				
				ktn[length]='\0';
				
			return ktn;
 }     
 
 
char * stringIneq(RAT * ieq, int noVars, int noIeq){
	int i , j, z;	
	double z2;
	char puffer[10];	
	char * c= "";		
	int column = noVars + 1;	
	for (i = 0; i < noIeq; i++){
		for (j = 0; j<noVars-1; j++) {     //bis x_{n-1}			
/*			z = ieq[(i)*column+j].num; 	sprintf(puffer, "%i", z);						
			c = konkatenation(c,puffer);    
			c = konkatenation(c,"/");				 				 
			
			z = ieq[(i)*column+j].den.i;			
			sprintf(puffer, "%i", z);			
			c = konkatenation(c,puffer);
			c = konkatenation(c,"*x");*/	
			#ifdef PRINT_DEBUG
			printf("Zaehler: %f\n",(double)ieq[(i)*column+j].num);
			printf("Nenner: %f\n",(double)ieq[(i)*column+j].den.i);
			#endif
			z2 = (ieq[(i)*column+j].num==0)?0:((double)ieq[(i)*column+j].num/(double)ieq[(i)*column+j].den.i);
			#ifdef PRINT_DEBUG
			printf("Ergebnisf: %f\n",z2);
			printf("Ergebnisg: %g\n\n",z2);
			#endif
						
			sprintf(puffer, "%g", z2);			
			c = konkatenation(c,puffer);
			c = konkatenation(c,"*x");
			
			z = j; 	sprintf(puffer, "%i", z); 				
			c = konkatenation(c,puffer);
			c = konkatenation(c," + ");
						
		}		
/*		z = ieq[(i)*column+noVars-1].num; 	 sprintf(puffer, "%d", z);				
		c = konkatenation(c,puffer);		          
		c = konkatenation(c,"/");			           
			
		z = ieq[(i)*column+noVars-1].den.i; 	sprintf(puffer, "%d", z); 			
		c = konkatenation(c,puffer);			             
		c = konkatenation(c,"*x");*/                           
		#ifdef PRINT_DEBUG
		printf("Zaehler: %f\n",(double)ieq[(i)*column+noVars-1].num);
		printf("Nenner: %f\n",(double)ieq[(i)*column+noVars-1].den.i);
		#endif
		z2 = (ieq[(i)*column+noVars-1].num==0)?0:((double)ieq[(i)*column+noVars-1].num/(double)ieq[(i)*column+noVars-1].den.i); 	 
		#ifdef PRINT_DEBUG
		printf("Ergebnisf: %f\n",z2);
		printf("Ergebnisg: %g\n\n",z2);
		#endif

		sprintf(puffer, "%g", z2);				
		c = konkatenation(c,puffer);		          
		c = konkatenation(c,"*x");                           

		z = noVars-1; 	sprintf(puffer, "%d", z); 				
		c = konkatenation(c,puffer);			                   
		c = konkatenation(c,"<=");
			
/*		z = ieq[(i)*column+noVars].num; 	sprintf(puffer, "%d", z); 			
		c = konkatenation(c,puffer);	
		c = konkatenation(c,"/");
				
		z = ieq[(i)*column+noVars].den.i; 	sprintf(puffer, "%d", z);		 						
		c = konkatenation(c,puffer);				
		c = konkatenation(c,"\n");*/				

		z2 = (ieq[(i)*column+noVars].num==0)?0:((double)ieq[(i)*column+noVars].num/(double)ieq[(i)*column+noVars].den.i); 	
	
		#ifdef PRINT_DEBUG
		printf("Zaehler: %f\n",(double)ieq[(i)*column+noVars].num);
		printf("Nenner: %f\n",(double)ieq[(i)*column+noVars].den.i);
		#endif
		sprintf(puffer, "%g", z2); 			
		#ifdef PRINT_DEBUG
		printf("Ergebnisf: %f\n",z2);
		printf("Ergebnisg: %g\n\n",z2);
		#endif
		c = konkatenation(c,puffer);				
		c = konkatenation(c,", ");				
		
	}
			//printf(c);
			return c;
}


 


	

void printMatrix(double * m, int rows,int columns){
	int i,j;
	for (i = 0; i < rows; i++){
		for (j = 0; j < columns; j++){
			printf(" %f ",m[i*columns+j]);
		}
		printf("\n");
	}
}


void printMatrix2v(double * m, int * rc){
	int i,j;
	for (i = 0; i < rc[0]; i++){
		for (j = 0; j < rc[1]; j++){
			printf(" %f ",m[i*rc[1]+j]);
		}
		printf("\n");
	}
}

