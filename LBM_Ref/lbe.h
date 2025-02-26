/*  Lattice-Boltzmann header file: 19 velocity model  */
/***********************************************************************
Susp3D: Lattice-Boltzmann simulation code for particle-fluid suspensions
Copyright (C) 2003 Tony Ladd

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
***********************************************************************/

#define TRUE       1
#define FALSE      0
#define DIMS       3
#define Rho_Fl    36
#define Num_Dir   19
#define Num_Dir_X  5
#define Num_Dir_NConserve 15
#define CS2       (1.0/3.0) 

static int x_p[ 5] = { 1,11,14,15,17};
static int x_m[ 5] = { 2,12,13,16,18};
static int y_p[ 5] = { 3, 7, 9,15,18};
static int y_m[ 5] = { 4, 8,10,16,17};
static int z_p[ 5] = { 5, 7,10,11,13};
static int z_m[ 5] = { 6, 8, 9,12,14};

/*  Directions         0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18  */
static int fac[19] = {12, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
static int c_x[19] = { 0, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1,-1,-1, 1, 1,-1, 1,-1};
static int c_y[19] = { 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 1,-1,-1, 1};
static int c_z[19] = { 0, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 1,-1, 1,-1, 0, 0, 0, 0};

/********************************************************
 19 modes for the unnormalized eigenvector e_q;
0-3 :   1, c_x, c_y, c_z
4-9 :   c^2-1, 3c_x^2 - c^2, c_y^2-c_z^2, c_x c_y, c_y c_z, c_z c_x
10-18 : (3c^2-5)c_x, (3c^2-5)c_y, (3c^2-5)c_z
        (c_y^2-c_z^2)c_x, (c_z^2-c_x^2)c_y, (c_x^2-c_y^2)c_z
        3c^4-6c^2+1, (2c^2-3)(3c_x^2-c^2), (2c^2-3)(c_y^2-c_z^2)


normalized eigenvector ne_qk = e_qk sqrt(fac[k]/w_k[q])
	**************************************************/

/* evector */
Float evector[Num_Dir][Num_Dir];
Float nmevector[Num_Dir][Num_Dir];
Float lin_evector[(Num_Dir+1)*Num_Dir];

/* normalization factor for eigenvector */
int wk[19];

/* inverse(evector) */ 
Float inv_evector[Num_Dir][Num_Dir];
Float inv_nmevector[Num_Dir][Num_Dir];
Float lin_invevector[(Num_Dir+1)*Num_Dir];


