/*
 * This file is part of RPorta. For licensing and copyright information
 * please see the file COPYING in the root directory of this
 * distribution or contact <robin.nunkesser@tu-dortmund.de>.
 * 
 * This file is a modification of the original file distributed with
 * PORTA (http://www.zib.de/Optimization/Software/Porta/).
 * Last modification: $Date: 2009/02/03 16:00:03 $
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
 

FILENAME: valid.h

AUTHOR: Thomas Christof

REVISED BY MECHTHILD STOER

REVISED BY ANDREAS LOEBEL
           ZIB BERLIN
           TAKUSTR.7
           D-14195 BERLIN

*******************************************************************************/
#define COMPILE_VALID
#ifdef COMPILE_VALID


#ifndef _VALID_H
#define _VALID_H


#include "porta.h"

int validmain( int argc,  char *argv[] );
extern int valid_points( int ieqNumber, int ineqOrEq, int, RAT *, int, int, RAT *, int, int, int, char ** );
extern void valid_ieqs( int, RAT *, int, int *, int *, int, RAT *, int, int, char ** );
extern void valid_ints( int, RAT *, int, int, int, RAT *, int, int, char * );
void init_construct_fcpt(int ineq);
void write_poi_file_fctp( int ineqNumber, int ineqOrEq, int dim, int lr, int flr, 
        int cone, int fce, int conv, int fcv );

#endif // _VALID_H
#endif
