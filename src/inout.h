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
 

FILENAME: inout.h

AUTHOR: Thomas Christof

REVISED BY MECHTHILD STOER

REVISED BY ANDREAS LOEBEL
           ZIB BERLIN
           TAKUSTR.7
           D-14195 BERLIN

*******************************************************************************/

#include <R.h>
#include <Rinternals.h>  
#include <Rdefines.h>
SEXP IEQPOI;
#ifndef _INOUT_H
#define _INOUT_H

 
#include "porta.h"


extern int scan_line2( int, char [], char *, char [] );
extern int read_input_file( char *, FILE *, int *, RAT **, int *, char *, int **, char *,
                            int **, char *, RAT ** );
extern void read_eqie( RAT **, int, int *, int *, int *, int *, char *, char [], char * );
extern void write_ieq_file( char *, FILE *, int, int, int, int *, 
                            int, int, int, int * );
extern void write_poi_file( char *, FILE *, int, int, int, int, int, int, int );
extern void writepoionie( FILE *, int, int, int, int );
extern void writesys( FILE *, int, int, int, int, int *, char, int * );
extern void writemark( FILE *, unsigned *, int, int * );
extern void writestatline( FILE *, int * );
extern FILE *wfopen( char * );
extern void I_RAT_writeline( FILE *, int, RAT *, int, RAT *, char, int * );
extern void constructIeqSexp(int dim,int n_ineqeq, RAT * ar1,int nel_ar1, RAT * ar6, int nel_ar6, int * elim_ord);
extern void constructPoiSexp(int dim,int n_points, RAT * ar1,int nel_ar1);

#endif // _INOUT_H
