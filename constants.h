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

#endif
