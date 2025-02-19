/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                            David Alonso, 2025 (c)                         */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <Include/MODEL.h>

#include "global.h"

gsl_rng * r; /* Global generator defined in main.c */

/* This code calculates the stochastic and determinisitic temporal evolution of several resource-consumer models 

   Compilation (see makefile environment variable 'MODEL'). Some compilation commands as example:

   . ~$ make
   . ~$ make STATIONARITY=STATIONARY_POINT_REPRESENTATION MODEL=DIFFUSION_HII_nD
   . ~$ make STATIONARITY=NON_STATIONARY_POINT_REPRESENTATION MODEL=DIFFUSION_HII_nD

   Exectution:
   
   MODEL = DIFFUSION_HII_nD
   
   Feeding experiments at a constant number of total consumers (-HN 20): 
   .~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 3 -HM 1 -HX 1 -HY 1 \
                          -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 \
                          -tn 50 -t0 0.0 -t1 1.5 -t4 0 -tR 10 -xn 0 -xN 20.0 \
                          -G2 1 -G3 0.0 -G4 1.5 -G5 1 -G6 0.0 -G7 20 \
                          -HK 10000 -HuR 0.0 -HuC 0.0 -H0 5.0 -H2 1.0 -H5 0.0 \
                          -H9 2.5 -H10 10.0 -Hp1 0.3725 -Hp2 0.5 -HN 20

   .~$ ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 3 -HM 1 -HX 1 -HY 1 \
                                 -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 \
                                 -tn 10 -t0 0.0 -t1 2.5 -t4 0 -tR 100 -xn 0 -xN 20.0 \
                                 -G2 1 -G3 0.0 -G4 2.5 -G5 1 -G6 0.0 -G7 8 \
                                 -HK 10000 -HuR 0.0 -HuC 0.0 -H0 5.0 -H2 1.0 -H5 0.0 \
                                 -H9 10.5 -H10 10.0 -Hp1 0.3725 -Hp2 1.0 -HN 20 -tE 0.2

   .~S ./DIFFUSION_HII_nD -y0 16 -y2 1 -HS 3 -HM 1 -HX 1 -HY 1 \
                          -n 3 -v0 0 -v1 1 -v2 2 -G0 1 -G1 3 \
                          -tn 10 -t0 0.0 -t1 2.5 -t4 0 -tR 100 -xn 0 -xN 20.0 \
                          -G2 1 -G3 0.0 -G4 2.5 -G5 1 -G6 0.0 -G7 20 \
                          -HK 10000 -HuR 0.0 -HuC 0.0 -H0 5.0 -H2 1.0 -H5 0.0 \
                          -H9 10.5 -H10 0.1  -Hp1 0.3725 -Hp2 1.0 -HN 20 -tE 0.2
   
Relevant input arguments for model DIFFUSION_HII_nD:
   -HK  [N or Total Carrying Capacity (for resources)]
   -Hp1 [approx y_R = f_i * p1 * N: Different resource levels at time 0.0, with f_i = 0.5, 0.3, 0.2, but this can be changed] 
   -Hp2 [Total No of Free Consumers at time 0.0: p2 * No_of_CONSUMERS]
   -HN 20 [Total No of CONSUMERS]
   -H9 and -H10 [Alpha_C_0 and Nu_C_0 for 1st Resource Type. Alpha_C_0 is the same for all resource types, but this can be changed.  
   -H0 and -H2 are Lambda_R_0 and Lambda_R1, which are overloaded to create extra handling rates (Nu), for the 2nd and 3rd resource type. 
   Notice that for this model -H5 should be always zero (No immigration of consumers: No of CONSUMERS is constant). 
   
In general:
   -HuR -HuC are the jumping rates (only relevant if there are more than one cell or patch in the system)
   -H0 -H2 -H5 are the external immigration (Lambda_R_0, Lambda_R_1 and Lambda_C_0)
   -H20 is the establishment rate 
   -H1 -H3 are the death rates, Delta_R_0, Delta_R_1, for propagules/resources. 
   -H6  is Delta_C_0, the death rate for both searching and handling consumers.
   -H9  and -H10 are the Alpha_C_0 and Nu_C_0  Holling Type II model parameters 
*/

int main(int argc, char **argv)
{
  int i;
  Parameter_Table Table;
  Time_Control Time;
  Time_Dependence_Control Time_Dependence;
  P_ARG = &Table;

#include "default.c"

  /* Command line arguments */
  if(argc>1) ArgumentControl(argc,argv);

#include <include.Output_Variables.default.aux.c>
  
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___A_L_L_O_C(   &Table );
  P_A_R_A_M_E_T_E_R___T_A_B_L_E___U_P_L_O_A_D( &Table, Index_Output_Variables );
  printf(" Parameter_Table structure has been correctly allocated and initiated\n");

  /* B E G I N : Reserving memmory for Parameter Space */
#include <include.Parameter_Space.default.aux.c>
  if( No_of_PARAMETERS == Table.TOTAL_No_of_MODEL_PARAMETERS ) {
    /* Full parameter space is in place. See also Model_Variables_Code.c */
    for(i=0; i<Table.TOTAL_No_of_MODEL_PARAMETERS; i++) Index[i] = Table.Index[i];
    No_of_PARAMETERS = Table.TOTAL_No_of_MODEL_PARAMETERS;
  }
  Parameter_Space * Space = (Parameter_Space *)calloc(1, sizeof(Parameter_Space));
  Parameter_Space_Alloc( Space, No_of_PARAMETERS, d);
  Parameter_Space_Initialization( Space, No_of_PARAMETERS, TOLERANCE, MAX_No_of_ITERATIONS,
    d, Index, Ranges, Acc);
  Table.S = Space;
  printf(" Parameter_Space structure has been correctly allocated and initiated\n");
  /*     E N D : ------------------------------------- */

#include <gsl_random_number_Setup.c>
  // #if defined VERBOSE
  /* BEGIN: Checking Random Number Generator Setup */
  for(i=0; i<10; i++){
    printf( "f(%d)=%g, ", i, gsl_rng_uniform(r) );
    printf( "f_GAUS(%d)=%g\n", i, gsl_ran_gaussian(r, 1.0) );
  }
  printf("\n"); Print_Press_Key(1,0,".");
  /*   END: Checking Random Number Generator Setup */
  // #endif

  if (TYPE_of_TIME_DEPENDENCE == 0) {
    printf(" Time_Control structure will be allocated: \n");
    printf(" %d output variables of length %d points will be allocated\n",
	   SUB_OUTPUT_VARIABLES, I_Time);
    T_I_M_E___C_O_N_T_R_O_L___A_L_L_O_C( &Time, &Table, I_Time);
    T_I_M_E___C_O_N_T_R_O_L___U_P_L_O_A_D( &Time, &Table, I_Time);
    printf(" Time_Control structure has been correctly allocated and set up\n");
  }
  else {
    #include <include.Time_Dependence_Control.default.aux.c>
    printf(" Time_Dependence_Control and Time_Control structures will be allocated: \n");
    printf(" %d output variables of length %d points will be allocated\n",
    SUB_OUTPUT_VARIABLES, I_Time);
    Time_Dependence_Control_Alloc(&Time, &Time_Dependence, &Table,
				  I_Time, TIME_DEPENDENT_PARAMETERS, No_of_COVARIATES);

    int No_of_EMPIRICAL_TIMES = I_Time;
    // Number of columns in the data files of time-dependent parameters
    Time_Dependence_Control_Upload(&Time, &Time_Dependence, &Table,
				   I_Time, No_of_EMPIRICAL_TIMES,
				   TIME_DEPENDENT_PARAMETERS, TYPE_of_TIME_DEPENDENCE,
				   TYPE_0_PARAMETERS, TYPE_1_PARAMETERS, TYPE_2_PARAMETERS,
				   No_of_COVARIATES,
				   dependent_parameter, forcing_pattern,
				   "File_of_Covariates.dat", Name_of_FILE[0] );
    printf(" Both Time_Control and Time_Dependence_Control structures have been\n");
    printf(" correctly allocated and set up\n");
  }

#if defined CPGPLOT_REPRESENTATION
  Table.CPG = A_C_T_I_V_A_T_E___C_P_G_P_L_O_T ( SUB_OUTPUT_VARIABLES, I_Time, 0, CPG_DRIVER_NAME);
  // Table.CPG_STO = A_C_T_I_V_A_T_E___2nd___C_P_G_P_L_O_T (1, SUB_OUTPUT_VARIABLES, I_Time, 0, "/TPNG");
  Table.CPG_STO = A_C_T_I_V_A_T_E___2nd___C_P_G_P_L_O_T (1, SUB_OUTPUT_VARIABLES, I_Time, 0, CPG_DRIVER_NAME);
  
  printf(" Two Parameterh_CPGPLOT plotting structures have been correctly allocated and initiated\n");
  printf(" These will open two windows (or two ploting devices of the same kind)\n");
  printf(" Table.CPG will store deterministic dynamic variables to plot\n");
  printf(" Table.CPG_STO will store stochastic dynamic variables to plot\n");
  printf(" As a consquence, deterministic and stochastic dynamics can be plotted\n");
  printf(" on the same device to compare (as it is done here, indicated by the first\n");
  printf(" input argument (0) of the A_Ch_T_I_V_A_T_E___2nd___C_P_G_P_L_O_T function).\n");
  printf(" Alternatively, two different devices (two different pdf files, for instance)\n");
  printf(" can be used, if required (1).\n");
#endif

  if(Table.TYPE_of_MODEL == 16) { // DIFFUSION_HII_nD
    	/* Also model where the TOTAL_No_of_CONSUMERS is a CONSTANT */
    	/* and they feed on multiple resources                      */
    	Common_Initial_Condition_Command_Line_Arguments_into_Table(&Table);
    	Resetting_Alpha_Nu_Vectors (&Table);
    	Resetting_Multiresource_Levels (&Table);  
    	Writing_Alpha_Nu_Theta_Vectors(&Table);  
  }
  else{
  	printf("This code can only use TYPE_of_MODEL = 16\n");
  	printf("which correspond to the multi-resource Holling Type II model\n");
  }

  Parameter_Values_into_Parameter_Table(&Table);
  M_O_D_E_L( &Table );

  // Some models (such as DIFFUSION_1R1C_2D and so on) do not have a stochastic
  // counter-part implemented yet!

  /* Stochastic Time Dynamics: A number of stochastic realizations will be produced */
  int No_of_RESOURCES_ACTUAL = Table.No_of_RESOURCES; 
  Parameter_Values_into_Parameter_Table(&Table);
  M_O_D_E_L___S_T_O( &Table );

  /* BEGIN : -------------------------------------------------------------------------
   */
  char boundary_File[80];
  sprintf(boundary_File, "boundary_Model_Parameter.c");
  write_Parameter_Table___RANGES___VALUES___LATEX ( "Latex_Parameter_Table.tex",
                                                    boundary_File,
                                                    &Table,
                                                    Space->P_MAX->data,
                                                    Space->P_min->data, Space->No_of_PARAMETERS );
  /*  END : ------------------------------------------------------------------------*/

  /* BEGIN : Freeing All Memmory * * * * * * * * * * * * * * */
#if defined CPGPLOT_REPRESENTATION
  #include <include.CPG.default.free.c>
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG, SUB_OUTPUT_VARIABLES );
  P_A_R_A_M_E_T_E_R___C_P_G_P_L_O_T___F_R_E_E( Table.CPG_STO, SUB_OUTPUT_VARIABLES );
  cpgclos();
#endif

#include <include.Parameter_Space.default.free.c>
  Parameter_Space_Free(Space, No_of_PARAMETERS); free( Space );

#include <include.Initial_Conditions.default.free.c>

#include <include.Output_Variables.default.free.c>

#include <include.Time_Dependence_Control.default.free.c>
  if (TYPE_of_TIME_DEPENDENCE == 0) T_I_M_E___C_O_N_T_R_O_L___F_R_E_E( &Time, &Table );
  else                        Time_Dependence_Control_Free( &Time_Dependence, &Table );

  P_A_R_A_M_E_T_E_R___T_A_B_L_E___F_R_E_E( &Table );
  /*  END : Freeing  All Memmory * * * * * * * * * * * * * * */

  printf("\nEnd of progam\n");
  return (0);
}
  
