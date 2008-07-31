/*
 * This file is part of RPorta. For licensing and copyright information
 * please see the file COPYING in the root directory of this
 * distribution or contact <robin.nunkesser@tu-dortmund.de>.
 * 
 * This file is a modification of the original file distributed with
 * PORTA (http://www.zib.de/Optimization/Software/Porta/).
 * Last modification: $Date: 2008/08/06 11:46:39 $
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


#ifndef _ARITH_H
#define _ARITH_H


#include "porta.h"


extern void I_RAT_assign( RAT *, RAT * );
extern void I_RAT_add( RAT, RAT, RAT * );
extern void I_RAT_sub( RAT, RAT, RAT * );
extern void I_RAT_mul( RAT, RAT, RAT * );
extern void I_RAT_row_prim( RAT *, RAT *, RAT *, int );
extern void gauss_calcnewrow( RAT *, RAT *, int, RAT *, int, int );
extern void vecpr( RAT *, RAT *, RAT *, int );
extern void row_add( RAT *, RAT *, RAT *, int );
extern int eqie_satisfied( RAT *, RAT *, int, int );
extern void scal_mul( RAT *, RAT *, RAT *, int );
extern int igcd( int, int );
extern int gcdrow( int *, int );


#endif // _ARITH_H
