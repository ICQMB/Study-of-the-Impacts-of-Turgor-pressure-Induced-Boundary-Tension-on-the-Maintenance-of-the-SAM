//externs.h
//==========
//Include guards
#ifndef _EXTERNS_H_INCLUDED_
#define _EXTERNS_H_INCLUDED_
//Declares external nonconstant global variables that will be
//defined in phys.h
extern int DIV_MECHANISM;
extern bool WUS_LEVEL;
extern double WUS_RAD_CONTRACTION_FACTOR;
extern double OOP_PROBABILITY;
extern double CK_RAD_CONTRACTION_FACTOR;
extern double MECH_DIV_PROB;
extern bool OUT_OF_PLANE_GROWTH;
extern bool UNIFORM_OOP_DIST;
extern bool CHEMICAL_GD;
extern bool BOUNDARY_PULL;
extern int BOUNDARY_PULL_TYPE;
extern double BOUNDARY_FORCE_MAGNITUDE;
extern bool L1_L2_FORCED_ANTICLINAL_DIV;
extern double THETA_ABC;
//Periclinal exempt as an idea
extern int TENSILE_CALC;
extern int NUM_STEPS_PER_FRAME;
extern int RECENT_DIV_NUM_FRAMES;
extern int VTK_PER_DATA_POINT;
extern int Weird_WUS;
#endif
