#ifndef CLASSIC6DOFKINE_H_
#define CLASSIC6DOFKINE_H_


enum Kine6dConfig {
	RA_AE_UW,
	RA_AE_DW,
	RA_BE_UW,
	RA_BE_DW,
	LA_AE_UW,
	LA_AE_DW,
	LA_BE_UW,
	LA_BE_DW,
};

typedef struct Kine6d_ Kine6d;
struct Kine6d_{
	float X, Y, Z;
	float A, B, C;
	float R[9];
	int fgR;
};

typedef struct Kine6dSol_ Kine6dSol;
struct Kine6dSol_ {
	float sol[8][6];
	unsigned char sol_flag[8][3];
};

void classic6dofForKine(float* q_, Kine6d* pose_);

void classic6dofInvKine(Kine6d* pose_, float* q_last_, Kine6dSol* q_);

#endif