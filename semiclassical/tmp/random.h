/*
 * Copyright (C) 2003 Klaus Frahm <frahm@irsamc.ups-tlse.fr>
 * Quantware MIPS Center, Laboratoire de Physique Theorique
 * University Paul Sabatier, Toulouse III
 * 118, route de Narbonne, 31062 Toulouse Cedex 4 - FRANCE
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 2 as
 * published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public Licens
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-
 *
 */

#ifndef _RANDOM

#define _RANDOM


void reset_random(void);
void save_seed(char *name=NULL);
int load_seed(char *name=NULL);
void new_random(void);
int rand_eig(void);
double rand_double(void);
double random_rectangular_old(double w);
double random_rectangular(double w);
double random_gauss(double w,int wn);

#endif /* RANDOM */
