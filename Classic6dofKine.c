#include "Classic6dofKine.h"
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <string.h>
#include <math.h>

#define sqrtf(X) (float)sqrt(X)
#define sinf(X) (float)sin(X)
#define cosf(X) (float)cos(X)
#define asinf(X) (float)asin(X)
#define acosf(X) (float)acos(X)
#define atanf(X) (float)atan(X)
#define atan2f(Y, X) (float)atan2(Y, X)


#define CLASSIC6DOF_L_BS 0.100f
#define CLASSIC6DOF_D_BS 0.050f
#define CLASSIC6DOF_L_SE 0.340f
#define CLASSIC6DOF_L_EW 0.360f
#define CLASSIC6DOF_D_EW 0.030f
#define CLASSIC6DOF_L_WT 0.100f

static float classic6dof_DH[6][4] = {
	{	0.0f,			CLASSIC6DOF_L_BS,	CLASSIC6DOF_D_BS,	-(float)M_PI_2	},
	{	-(float)M_PI_2,	0.0f,				CLASSIC6DOF_L_SE,	0.0f			},
	{	 (float)M_PI_2,	CLASSIC6DOF_D_EW,	0.0f,				 (float)M_PI_2	},
	{	0.0f,			CLASSIC6DOF_L_EW,	0.0f,				-(float)M_PI_2	},
	{	0.0f,			0.0f,				0.0f,				 (float)M_PI_2	},
	{	0.0f,			CLASSIC6DOF_L_WT,	0.0f,				0.0f			}
}; // home, d, a, alpha

static float L1_bs[3] =	{	 CLASSIC6DOF_D_BS,	-CLASSIC6DOF_L_BS,	0.0f				};
static float L2_se[3] = {	 CLASSIC6DOF_L_SE,	 0.0f,				0.0f				};
static float L3_ew[3] = { -CLASSIC6DOF_D_EW,	 0.0f,				CLASSIC6DOF_L_EW	};
static float L6_wt[3] = {	 0.0f,				 0.0f,				CLASSIC6DOF_L_WT	};


static void matMultiply(float* M1, float* M2, float* M, int m, int l, int n)
{
	float tmp;
	int i, j, k;
	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			tmp = 0.0f;
			for (k = 0; k < l; k++) {
				tmp += M1[l*i + k] * M2[n*k+j];
			}
			M[n*i+j] = tmp;
		}
	}
}

static void matRotMatToFixedAngle(float* R, float* fa)
{
	float A, B, C, cb;
	if (fabs(R[6]) >= 1.0 - 0.0001) {
		if (R[6] < 0) {
			A = 0.0f;
			B =  (float)M_PI_2;
			C = atan2f(R[1], R[4]);
		} else {
			A = 0.0f;
			B = -(float)M_PI_2;
			C = -atan2f(R[1], R[4]);
		}
	} else {
		B = atan2f(-R[6], sqrtf(R[0] * R[0] + R[3] * R[3]));
		cb = cosf(B);
		A = atan2f(R[3] / cb, R[0] / cb);
		C = atan2f(R[7] / cb, R[8] / cb);
	}
	fa[0] = C;
	fa[1] = B;
	fa[2] = A;
}
static void matFixedAngleToRotMat(float* fa, float* R)
{
	float ca, cb, cc, sa, sb, sc;
	cc = cosf(fa[0]);
	cb = cosf(fa[1]);
	ca = cosf(fa[2]);
	sc = sinf(fa[0]);
	sb = sinf(fa[1]);
	sa = sinf(fa[2]);

	R[0] = ca*cb; R[1] = ca*sb*sc - sa*cc; R[2] = ca*sb*cc + sa*sc;
	R[3] = sa*cb; R[4] = sa*sb*sc + ca*cc; R[5] = sa*sb*cc - ca*sc;
	R[6] = -sb;   R[7] = cb*sc;            R[8] = cb*cc;
}

void classic6dofForKine(float* q_, Kine6d* pose_)
{
	float q[6];
	float cosq, sinq;
	float cosa, sina;
	float d, a;
	float P06[6];
	float R06[9];
	float R[6][9];
	float R02[9];
	float R03[9];
	float R04[9];
	float R05[9];
	float L0_bs[3];
	float L0_se[3];
	float L0_ew[3];
	float L0_wt[3];
	
	int i;

	for (i = 0; i < 6; i++) {
		q[i] = q_[i] + classic6dof_DH[i][0];
		cosq = cosf(q[i]);
		sinq = sinf(q[i]);
		cosa = cosf(classic6dof_DH[i][3]);
		sina = sinf(classic6dof_DH[i][3]);
		d = classic6dof_DH[i][1];
		a = classic6dof_DH[i][2];

		R[i][0] = cosq; R[i][1] = -cosa * sinq; R[i][2] =  sina * sinq;
		R[i][3] = sinq; R[i][4] =  cosa * cosq; R[i][5] = -sina * cosq;
		R[i][6] = 0.0f; R[i][7] =  sina;        R[i][8] =  cosa;
	}

	matMultiply(R[0], R[1], R02, 3, 3, 3);
	matMultiply(R02 , R[2], R03, 3, 3, 3);
	matMultiply(R03 , R[3], R04, 3, 3, 3);
	matMultiply(R04 , R[4], R05, 3, 3, 3);
	matMultiply(R05 , R[5], R06, 3, 3, 3);

	matMultiply(R[0], L1_bs, L0_bs, 3, 3, 1);
	matMultiply(R02 , L2_se, L0_se, 3, 3, 1);
	matMultiply(R03 , L3_ew, L0_ew, 3, 3, 1);
	matMultiply(R06 , L6_wt, L0_wt, 3, 3, 1);

	for (i = 0; i < 3; i++) {
		P06[i] = L0_bs[i] + L0_se[i] + L0_ew[i] + L0_wt[i];
	}
	matRotMatToFixedAngle(R06, &P06[3]);

	pose_->X = P06[0];
	pose_->Y = P06[1];
	pose_->Z = P06[2];
	pose_->A = P06[3];
	pose_->B = P06[4];
	pose_->C = P06[5];
	memcpy(pose_->R, R06, 9*sizeof(float));
}

void classic6dofInvKine(Kine6d* pose_, float* q_last_, Kine6dSol* q_)
{
	static float l_se_2 = CLASSIC6DOF_L_SE * CLASSIC6DOF_L_SE;
	static float l_se = CLASSIC6DOF_L_SE;
	static float l_ew_2 = CLASSIC6DOF_L_EW * CLASSIC6DOF_L_EW + CLASSIC6DOF_D_EW * CLASSIC6DOF_D_EW;
	static float l_ew = 0;
	static float atan_e = 0;

	float qs[2];
	float qa[2][2];
	float qw[2][3];
	float cosqs, sinqs;
	float cosqa[2], sinqa[2];
	float cosqw, sinqw;
	float P06[6];
	float R06[9];
	float P0_w[3];
	float P1_w[3];
	float L0_wt[3];
	float L1_sw[3];
	float R10[9];
	float R31[9];
	float R30[9];
	float R36[9];
	float l_sw_2, l_sw, atan_a, acos_a, acos_e;

	int ind_arm, ind_elbow, ind_wrist;
	int i;

	if (0 == l_ew) {
		l_ew = sqrtf(l_ew_2);
		atan_e = atanf(CLASSIC6DOF_D_EW / CLASSIC6DOF_L_EW);
	}

	P06[0] = pose_->X;
	P06[1] = pose_->Y;
	P06[2] = pose_->Z;
	if (0 == pose_->fgR) {
		P06[3] = pose_->A;
		P06[4] = pose_->B;
		P06[5] = pose_->C;
		matFixedAngleToRotMat(&P06[3], R06);
	} else {
		memcpy(R06, pose_->R, 9*sizeof(float));
	}
	for (i = 0; i < 2; i++) {
		qs[i] = q_last_[0];
		qa[i][0] = q_last_[1]; qa[i][1] = q_last_[2];
		qw[i][0] = q_last_[3]; qw[i][1] = q_last_[4]; qw[i][2] = q_last_[5];
	}
	// q1 solution pair ///////////
	matMultiply(R06, L6_wt, L0_wt, 3, 3, 1);
	for (i = 0; i < 3; i++) {
		P0_w[i] = P06[i] - L0_wt[i];
	}
	if (sqrt(P0_w[0]*P0_w[0] + P0_w[1]*P0_w[1]) <= 0.000001) {
		qs[0] = q_last_[0]; // right arm
		qs[1] = q_last_[0]; // left arm
		for (i = 0; i < 4; i++) {
			q_->sol_flag[0 + i][0] = -1;
			q_->sol_flag[4 + i][0] = -1;
		}
	} else {
		qs[0] = atan2f( P0_w[1],  P0_w[0]); // right arm
		qs[1] = atan2f(-P0_w[1], -P0_w[0]); // left arm
		for (i = 0; i < 4; i++) {
			q_->sol_flag[0 + i][0] =  1;
			q_->sol_flag[4 + i][0] =  1;
		}
	}
	// two arm config. ////////////
	for (ind_arm = 0; ind_arm < 2; ind_arm++) {
		// q2, q3 solution pair ///
		cosqs = cosf(qs[ind_arm] + classic6dof_DH[0][0]);
		sinqs = sinf(qs[ind_arm] + classic6dof_DH[0][0]);

		R10[0] =  cosqs; R10[1] = sinqs; R10[2] =  0.0f;
		R10[3] =   0.0f; R10[4] =  0.0f; R10[5] = -1.0f;
		R10[6] = -sinqs; R10[7] = cosqs; R10[8] =  0.0f;

		matMultiply(R10, P0_w, P1_w, 3, 3, 1);
		for (i = 0; i < 3; i++) {
			L1_sw[i] = P1_w[i] - L1_bs[i];
		}
		l_sw_2 = L1_sw[0]*L1_sw[0] + L1_sw[1]*L1_sw[1];
		l_sw = sqrtf(l_sw_2);

		if			(fabs(l_se + l_ew - l_sw) <= 0.000001) {
			qa[0][0] = atan2f(L1_sw[1], L1_sw[0]);
			qa[1][0] = qa[0][0];
			qa[0][1] = 0.0f;
			qa[1][1] = 0.0f;
			if (l_sw > l_se + l_ew) {
				for (i = 0; i < 2; i++) {
					q_->sol_flag[4*ind_arm + 0 + i][1] = 0;
					q_->sol_flag[4*ind_arm + 2 + i][1] = 0;
				}
			} else {
				for (i = 0; i < 2; i++) {
					q_->sol_flag[4*ind_arm + 0 + i][1] = 1;
					q_->sol_flag[4*ind_arm + 2 + i][1] = 1;
				}
			}
		} else if	(fabs(l_sw - fabs(l_se - l_ew)) <= 0.000001) {
			qa[0][0] = atan2f(L1_sw[1], L1_sw[0]);
			qa[1][0] = qa[0][0];
			if (0 == ind_arm) { // right arm
				qa[0][1] =  (float)M_PI; // above elbow
				qa[1][1] = -(float)M_PI; // below elbow
			} else { // /////// // left arm
				qa[0][1] = -(float)M_PI; // above elbow
				qa[1][1] =  (float)M_PI; // below elbow
			}
			if	(l_sw < fabs(l_se - l_ew)) {
				for (i = 0; i < 2; i++) {
					q_->sol_flag[4*ind_arm + 0 + i][1] = 0;
					q_->sol_flag[4*ind_arm + 2 + i][1] = 0;
				}
			} else {
				for (i = 0; i < 2; i++) {
					q_->sol_flag[4*ind_arm + 0 + i][1] = 1;
					q_->sol_flag[4*ind_arm + 2 + i][1] = 1;
				}
			}
		} else {
			atan_a = atan2f(L1_sw[1], L1_sw[0]);
			acos_a = 0.5f*(l_se_2 + l_sw_2 - l_ew_2) / (l_se*l_sw);
			if	(acos_a >=  1.0f) acos_a = 0.0f;
			else if	(acos_a <= -1.0f) acos_a = (float)M_PI;
			else	acos_a = acosf(acos_a);
			acos_e = 0.5f*(l_se_2 + l_ew_2 - l_sw_2) / (l_se*l_ew);
			if	(acos_e >=  1.0f) acos_e = 0.0f;
			else if	(acos_e <= -1.0f) acos_e = (float)M_PI;
			else	acos_e = acosf(acos_e);
			if (0 == ind_arm) { // right arm
				// above elbow
				qa[0][0] = atan_a - acos_a + (float)M_PI_2;
				qa[0][1] = atan_e - acos_e + (float)M_PI;
				// below elbow
				qa[1][0] = atan_a + acos_a + (float)M_PI_2;
				qa[1][1] = atan_e + acos_e - (float)M_PI;

			} else { // /////// // left arm
				// above elbow
				qa[0][0] = atan_a + acos_a + (float)M_PI_2;
				qa[0][1] = atan_e + acos_e - (float)M_PI;
				// below elbow
				qa[1][0] = atan_a - acos_a + (float)M_PI_2;
				qa[1][1] = atan_e - acos_e + (float)M_PI;
			}
			for (i = 0; i < 2; i++) {
				q_->sol_flag[4*ind_arm + 0 + i][1] = 1;
				q_->sol_flag[4*ind_arm + 2 + i][1] = 1;
			}
		}
		// two elbow config. ////////
		for (ind_elbow = 0; ind_elbow < 2; ind_elbow++) {
			// q3,q4,q5 solution pair
			cosqa[0] = cosf(qa[ind_elbow][0] + classic6dof_DH[1][0]); sinqa[0] = sinf(qa[ind_elbow][0] + classic6dof_DH[1][0]);
			cosqa[1] = cosf(qa[ind_elbow][1] + classic6dof_DH[2][0]); sinqa[1] = sinf(qa[ind_elbow][1] + classic6dof_DH[2][0]);

			R31[0] = cosqa[0]*cosqa[1] - sinqa[0]*sinqa[1]; R31[1] =   cosqa[0]*sinqa[1] + sinqa[0]*cosqa[1]; R31[2] = 0.0f;
			R31[3] = 0.0f; R31[4] = 0.0f; R31[5] = 1.0f;
			R31[6] = cosqa[0]*sinqa[1] + sinqa[0]*cosqa[1]; R31[7] = - cosqa[0]*cosqa[1] + sinqa[0]*sinqa[1]; R31[8] = 0.0f;

			matMultiply(R31, R10, R30, 3, 3, 3);
			matMultiply(R30, R06, R36, 3, 3, 3);

			if			(R36[8] >= 1.0 - 0.000001) {
				cosqw =  1.0f;
				qw[0][1] = 0.0f;
				qw[1][1] = 0.0f;
			} else if	(R36[8] <= -1.0 + 0.000001) {
				cosqw = -1.0f;
				if (0 == ind_arm) { // right arm
					qw[0][1] =  (float)M_PI;
					qw[1][1] = -(float)M_PI;
				} else { // /////// // left arm
					qw[0][1] = -(float)M_PI;
					qw[1][1] =  (float)M_PI;
				}
			} else {
				cosqw = R36[8];
				if (0 == ind_arm) { // right arm
					qw[0][1] =  acosf(cosqw); // up wrist
					qw[1][1] = -acosf(cosqw); // down wrist
				} else { // /////// // left arm
					qw[0][1] = -acosf(cosqw); // up wrist
					qw[1][1] =  acosf(cosqw); // down wrist
				}
			}
			if (1.0f == cosqw || -1.0f == cosqw) {
				if (0 == ind_arm) { // right arm
					// q4 = q_last
					qw[0][0] = q_last_[3];
					cosqw = cosf(q_last_[3] + classic6dof_DH[3][0]); sinqw = sinf(q_last_[3] + classic6dof_DH[3][0]);
					qw[0][2] = atan2f(cosqw*R36[3] - sinqw*R36[0], cosqw*R36[0] + sinqw*R36[3]);
					// q6 = q_last
					qw[1][2] = q_last_[5];
					cosqw = cosf(q_last_[5] + classic6dof_DH[5][0]); sinqw = sinf(q_last_[5] + classic6dof_DH[5][0]);
					qw[1][0] = atan2f(cosqw*R36[3] - sinqw*R36[0], cosqw*R36[0] + sinqw*R36[3]);
				} else { // /////// // left arm
					// q6 = q_last
					qw[0][2] = q_last_[5];
					cosqw = cosf(q_last_[5] + classic6dof_DH[5][0]); sinqw = sinf(q_last_[5] + classic6dof_DH[5][0]);
					qw[0][0] = atan2f(cosqw*R36[3] - sinqw*R36[0], cosqw*R36[0] + sinqw*R36[3]);
					// q4 = q_last
					qw[1][0] = q_last_[3];
					cosqw = cosf(q_last_[3] + classic6dof_DH[3][0]); sinqw = sinf(q_last_[3] + classic6dof_DH[3][0]);
					qw[1][2] = atan2f(cosqw*R36[3] - sinqw*R36[0], cosqw*R36[0] + sinqw*R36[3]);
				}
				q_->sol_flag[4*ind_arm+2*ind_elbow+0][2] = -1;
				q_->sol_flag[4*ind_arm+2*ind_elbow+1][2] = -1;
			} else {
				if (0 == ind_arm) { // right arm
					// q4
					qw[0][0] = atan2f( R36[5],  R36[2]); // up wrist
					qw[1][0] = atan2f(-R36[5], -R36[2]); // down wrist
					// q6
					qw[0][2] = atan2f( R36[7], -R36[6]); // up wrist
					qw[1][2] = atan2f(-R36[7],  R36[6]); // down wrist
				} else { // /////// // left arm
					// q4
					qw[0][0] = atan2f(-R36[5], -R36[2]); // up wrist
					qw[1][0] = atan2f( R36[5],  R36[2]); // down wrist
					// q6
					qw[0][2] = atan2f(-R36[7],  R36[6]); // up wrist
					qw[1][2] = atan2f( R36[7], -R36[6]); // down wrist
				}
				q_->sol_flag[4*ind_arm+2*ind_elbow+0][2] =  1;
				q_->sol_flag[4*ind_arm+2*ind_elbow+1][2] =  1;
			}
			// two wrist config. ////
			for (ind_wrist = 0; ind_wrist < 2; ind_wrist++) {
				if		(qs[ind_arm] >  (float)M_PI)
					q_->sol[4*ind_arm+2*ind_elbow+ind_wrist][0] = qs[ind_arm] - (float)M_PI;
				else if	(qs[ind_arm] < -(float)M_PI)
					q_->sol[4*ind_arm+2*ind_elbow+ind_wrist][0] = qs[ind_arm] + (float)M_PI;
				else
					q_->sol[4*ind_arm+2*ind_elbow+ind_wrist][0] = qs[ind_arm];
				for (i = 0; i < 2; i++) {
					if		(qa[ind_elbow][i] >  (float)M_PI)
						q_->sol[4*ind_arm+2*ind_elbow+ind_wrist][1 + i] = qa[ind_elbow][i] - (float)M_PI;
					else if	(qa[ind_elbow][i] < -(float)M_PI)
						q_->sol[4*ind_arm+2*ind_elbow+ind_wrist][1 + i] = qa[ind_elbow][i] + (float)M_PI;
					else
						q_->sol[4*ind_arm+2*ind_elbow+ind_wrist][1 + i] = qa[ind_elbow][i];
				}
				for (i = 0; i < 3; i++) {
					if		(qw[ind_wrist][i] >  (float)M_PI) 
						q_->sol[4*ind_arm+2*ind_elbow+ind_wrist][3 + i] = qw[ind_wrist][i] - (float)M_PI;
					else if	(qw[ind_wrist][i] < -(float)M_PI)
						q_->sol[4*ind_arm+2*ind_elbow+ind_wrist][3 + i] = qw[ind_wrist][i] + (float)M_PI;
					else
						q_->sol[4*ind_arm+2*ind_elbow+ind_wrist][3 + i] = qw[ind_wrist][i];
				}
			} // for ind_wrist
		} // for ind_elbow
	} // for ind_arm
	
}
