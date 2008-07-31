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
 

FILENAME: mp.h

AUTHOR: Thomas Christof

REVISED BY MECHTHILD STOER

REVISED BY ANDREAS LOEBEL
           ZIB BERLIN
           TAKUSTR.7
           D-14195 BERLIN

*******************************************************************************/


#ifndef _MP_H
#define _MP_H


#include "porta.h"


extern void arith_overflow_func( int, void (*)(), RAT, RAT, RAT * );
extern void size_info( RAT *, int *, int * );
extern void L_RAT_assign( RAT *, RAT * );
extern void L_RAT_to_RAT( RAT *, int );
extern int vals_lt_MAXINT( RAT *, int );
extern void RAT_to_L_RAT( RAT *, int );
extern void L_RAT_add( RAT, RAT, RAT * );
extern void L_RAT_sub( RAT, RAT, RAT * );
extern void L_RAT_mul( RAT, RAT, RAT * );
extern void L_RAT_kue( loint *, loint * );
extern void L_RAT_row_prim( RAT *, RAT *, RAT *, int );
extern void L_RAT_writeline( FILE *, int, RAT *, int, char, int * );
extern loint lgcdrow( loint *, int );
extern void hexprint( FILE *, loint );
extern int return_from_mp( );


#endif // _MP_H
