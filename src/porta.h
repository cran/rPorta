/*
 * This file is part of RPorta. For licensing and copyright information
 * please see the file COPYING in the root directory of this
 * distribution or contact <robin.nunkesser@tu-dortmund.de>.
 * 
 * This file is a modification of the original file distributed with
 * PORTA (http://www.zib.de/Optimization/Software/Porta/).
 * Last modification: $Date: 2008/05/06 14:16:21 $
 */

/*******************************************************************************

Copyright (C) 1997-2002 Thomas Christof and Andreas Loebel
 
This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.
 
This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 59 Temple
Place, Suite 330, Boston, MA 02111-1307 USA
 

FILENAME: arith.h

AUTHOR: Thomas Christof

REVISED BY MECHTHILD STOER

REVISED BY ANDREAS LOEBEL
           ZIB BERLIN
           TAKUSTR.7
           D-14195 BERLIN

*******************************************************************************/


#ifndef _PORTA_H
#define _PORTA_H


/* 21.01.1994 include version information */
#define VERSION "1.4.0-20021125"


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
/* include <strings.h> */

/* unsigned is supposed to be 32 bits long */

struct RAT {
  long num;
  union INT {
    int i;
    unsigned *p;
  } den;
};
typedef struct RAT RAT;

int mp_state;   
#define MP_not_ready (mp_state == -1)
#define MP_ready (mp_state == 0)
#define MP_realised (mp_state == 1)
#define SET_MP_not_ready mp_state = -1
#define SET_MP_ready mp_state = 0
#define SET_MP_realised mp_state = 1

void 
  (*RAT_add)(),
  (*RAT_sub)(),
  (*RAT_mul)(),
  (*RAT_row_prim)(),
  (*RAT_assign)(),
  (*writeline)();

#define MAX_LEN_LINT 20

struct loint {
  int len,neg;
  unsigned val[MAX_LEN_LINT];
};
typedef struct loint loint;

struct lorat {
  loint num,den;
};
typedef struct lorat lorat;


#define FIRST_SYS_EL 50000 
#define INCR_SYS_EL 100000 
#define INCR_LIST 5000  
#define INCR_INSYS_ROW 1000 

#define CP (char *)
#define U (unsigned )
#define UP (unsigned *)
#define RP (RAT *)
#define szU (sizeof(unsigned))

struct list {
  RAT *sys; 
  unsigned *mark;
  RAT *ptr;
  }  ;


typedef struct list *listp;

extern char andreas[];
extern listp *porta_list;
     
RAT RAT_const[2],var[4];

 


RAT *ar1,*ar2,*ar3,*ar4,*ar5,*ar6;
long nel_ar1,nel_ar2,nel_ar3,nel_ar4,nel_ar5,nel_ar6;

int  maxlist, total_size;

int  dim,
     equa,    /* number of equalities */
     ineq,    /* number of inequalities */
     conv, 
     cone, 
     points,
     blocks;


FILE *fp,*prt;
char * RATallo();
char * allo();

/*  options  */

int option, allowed_options;
#define is_set(x) (option & x)
#define Protocol_to_file 1
#define Redundance_check 4
#define Validity_table_out 8
#define Statistic_of_coefficients 16
#define Chernikov_rule_off 32
#define Fmel 64
#define Dim 128
#define Sort 256
#define Cfctp 512
#define Posie 1024
#define Iespo 2048
#define Vint 4096
#define Traf 8192
#define Opt_elim 16384
#define Long_arithmetic 32768
#define Datei_Parsing 65536  // added by sebsastian ruthe soll die Option Parsen darstellen

RAT * testPorta(RAT * ppoints, int punktAnz, int punktDim,int * pnoOfIneq);
int portamain( int argc, char *argv[] );
extern int no_denom( int, int, int, int );
extern void reallocate( int, RAT ** );
extern void polarformat( RAT *, int *, int, RAT * );
extern int *check_and_reorder_elim_ord( int *, int * );
extern void resubst( RAT *, int, int [] );
extern void blow_up( RAT *, int, int *, int, int );
extern void backwsubst( RAT *, int, int );
extern void origin_add( int, RAT * );
extern void gentableau( RAT *, int, int *, int ** );
extern void reorder_var( int, RAT *, RAT **, int *, int *, int **, int ** );

extern float time_used( );
extern void init_total_time( );
extern float total_time( );


#endif // _PORTA_H
