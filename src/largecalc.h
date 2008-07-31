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
 

FILENAME: largecalc.h

AUTHOR: Thomas Christof

REVISED BY MECHTHILD STOER

REVISED BY ANDREAS LOEBEL
           ZIB BERLIN
           TAKUSTR.7
           D-14195 BERLIN

*******************************************************************************/


#ifndef _LARGECALC_H
#define _LARGECALC_H


#include "porta.h"


extern void lsub( loint, loint, loint * );
extern void lsubber( unsigned *, unsigned *, unsigned *, int, int, int * );
extern void ladd( loint, loint, loint * );
extern void ladder( unsigned *, unsigned *, unsigned *, int, int, int * );
extern void lmul( loint, loint, loint * );
extern void lmuller( unsigned *, unsigned *, unsigned *, int, int, int * );
extern void porta_ldiv( loint, loint, loint *, loint * );
extern int lorder( unsigned *, unsigned *, int, int );
extern void lgcd( loint, loint, loint * );


#endif // _LARGECALC_H
