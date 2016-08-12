#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "Classic6dofKine.h"

int main()
{
	float q0[6] = { 0,0,0,0,0,0 };
	float q[6] = { 0,0,1.570796f,0,1.570796f,0 };
	float q_last[6] = { 0,0,0,0,0,0 };
	Kine6d pose;
	Kine6dSol q_sol;
	char ch = '\0';
	int i;

	pose.fgR = 0;

	classic6dofForKine(q0, &pose);
	printf("FK(q0)  :=< %.6f, %.6f, %.6f, %.6f, %.6f, %.6f >\n", pose.X, pose.Y, pose.Z, pose.A, pose.B, pose.C);

	while (1) {
		printf("quit 'q' or continue 'c' ?");
		//fflush(stdin);
		scanf(" %c", &ch);
		printf("%c\n", ch);
		if ('q' == ch) break;
		while ((ch = getchar()) != '\n' && ch != EOF);

		printf("input q:=");
		scanf("%f %f %f %f %f %f", &q[0], &q[1], &q[2], &q[3], &q[4], &q[5]);
		while ((ch = getchar()) != '\n' && ch != EOF);
		
		for (i = 0; i < 6; i++) {
			q[i] *= (3.14159265f / 180.0f);
		}
		printf("\nq       :=< %.6f, %.6f, %.6f, %.6f, %.6f, %.6f >\n", q[0], q[1], q[2], q[3], q[4], q[5]);
		
		classic6dofForKine(q, &pose);
		printf("FK      :=< %.6f, %.6f, %.6f, %.6f, %.6f, %.6f >\n", pose.X, pose.Y, pose.Z, pose.A, pose.B, pose.C);
		
		classic6dofInvKine(&pose, q_last, &q_sol);

		for (i = 0; i < 8; i++) {
			printf("q[%d]    :=< %.6f, %.6f, %.6f, %.6f, %.6f, %.6f >\n", i, q_sol.sol[i][0], q_sol.sol[i][1], q_sol.sol[i][2], q_sol.sol[i][3], q_sol.sol[i][4], q_sol.sol[i][5]);
			classic6dofForKine(q_sol.sol[i], &pose);
			printf("FK(q[%d]):=< %.6f, %.6f, %.6f, %.6f, %.6f, %.6f >\n", i, pose.X, pose.Y, pose.Z, pose.A, pose.B, pose.C);
		}

		memcpy(q_last, q, 6*sizeof(float));
	}
	//system("pause");
	return 0;
}