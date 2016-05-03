#ifndef CONSTANTS_H
#define CONSTANTS_H

//peak bin intensity weighting
#define PRIMARY_INTENSITY 50
#define FLANKING_INTENSITY 25

// constants
#define PROTON_MASS     1.00727646688
/* #define H_MASS          1.007825032 */
/* #define C_MASS          12.0 */
/* #define N_MASS          14.003074007 */
/* #define O_MASS          15.994914622 */

/* #define OH_MASS         17.002739654 */
/* #define NH3_MASS        17.026549103 */
/* #define H2O_MASS        18.010564686 */
/* #define CO_MASS         27.994914622 */
/* #define H3PO4_MASS      97.97689509 */

#define H_MASS 1.007825035
#define C_MASS 12.00000000
#define O_MASS 15.99491463
#define N_MASS 14.003074

#define OH_MASS  (O_MASS + H_MASS)
#define NH3_MASS (N_MASS + H_MASS*3)
#define H2O_MASS (H_MASS*2 + O_MASS)
#define CO_MASS  (C_MASS + O_MASS)

// amino acid monoisotopic masses
#define A_AA_MASS  71.03711378
#define B_AA_MASS 114.53493523 //Asn|Asp
#define C_AA_MASS 103.00918451
#define D_AA_MASS 115.02694302
#define E_AA_MASS 129.04259308
#define F_AA_MASS 147.06841390
#define G_AA_MASS  57.02146372
#define H_AA_MASS 137.05891186
#define I_AA_MASS 113.08406396
#define J_AA_MASS   0.0       
#define K_AA_MASS 128.09496300
#define L_AA_MASS 113.08406396
#define M_AA_MASS 131.04048463
#define N_AA_MASS 114.04292744
#define O_AA_MASS 237.147727   //Pyrrolysine
#define P_AA_MASS  97.05276384
#define Q_AA_MASS 128.05857750
#define R_AA_MASS 156.10111102
#define S_AA_MASS  87.03202840
#define T_AA_MASS 101.04767846
#define U_AA_MASS 150.953636   //Selenocysteine
#define V_AA_MASS  99.06841390
#define W_AA_MASS 186.07931294
#define X_AA_MASS   0.0       
#define Y_AA_MASS 163.06332852
#define Z_AA_MASS 128.55058529 //Glu|Gln

#endif
