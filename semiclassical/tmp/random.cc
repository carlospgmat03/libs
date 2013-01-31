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

#include <stdio.h>
#include "random.h"

#define root3 1.73205080756887719318
#define max_random 2147483647
#define max_zufall 2147483647
#define SEED_FILE ".random.seed"

/* Special random number generator using two seed numbers */

unsigned int konst_wert_1=217828199;
unsigned int konst_wert_2=314159269;
unsigned int seed_wert_1=111111;
unsigned int seed_wert_2=1111111;
unsigned int seed_old_1,seed_old_2;

void reset_random(){
  seed_wert_1=111111;
  seed_wert_2=1111111;
}

void save_seed(char *name){
  FILE *fp;

  if(name==NULL){
    fp=fopen(SEED_FILE,"w");
  }else{
    fp=fopen(name,"w");
  }
  fprintf(fp,"%12u %12u\n",seed_wert_1,seed_wert_2);
  fclose(fp);
}

int load_seed(char *name){
  FILE *fp;

  if(name==NULL){
    fp=fopen(SEED_FILE,"r");
  }else{
    fp=fopen(name,"r");
  }
  if(fp!=NULL){
    fscanf(fp,"%u%u",&seed_wert_1,&seed_wert_2);
    fclose(fp);
    return 1;
  }else{
    reset_random();
    return 0;
  }
}

void new_random(){
    printf("seed_1,  seed_2 ?  ");
    scanf("%u%u",&seed_wert_1,&seed_wert_2);
    printf("seed_1 = %12u ,  seed_2 = %12u\n",seed_wert_1,seed_wert_2);
    seed_old_1=seed_wert_1;
    seed_old_2=seed_wert_2;
}

int rand_eig(){
  unsigned int x;
  x=(konst_wert_1*seed_wert_2-konst_wert_2*seed_wert_1)%max_random;
  seed_wert_1=seed_wert_2;
  seed_wert_2=x;
  return x;
}
/* provides a random number between 0 and 1 */

double rand_double(void){
  int a;

  a=rand_eig();
  return (double)a/(double)max_zufall;
}

/* rectangular distribution between -w/2 and +w/2 */

double random_rectangular_old(double w){
  int a;

  a=rand_eig();
  return w*((double)a/(double)max_random-0.5);
}

/* rectangular distribution between -w and +w */

double random_rectangular(double w){
  int a;

  a=rand_eig();
/*  a=rand(); */
  return (w+w)*root3*((double)a/(double)max_random-0.5);
}

/* approximate gauss distribution of variance w */

double random_gauss(double w,int wn){
  double x;
  int i,n;

  x=0; n=wn*wn;
  for(i=0;i<n;i++) x+=random_rectangular(w);
  return x/(double)wn;
}

