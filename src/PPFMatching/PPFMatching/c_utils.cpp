
#include "c_utils.h"
#include <stdio.h>

void matrix_sum(int Am, int An, int Bm, int Bn, double *A,  double *B, double *R) 
{
	int r = Am;
	int c = An;
	int n = r * c, i;

	if (Am != Bm || An != Bn) 
	{
		printf("[matrix_sum] Error: mismatched dimensions\n");
		return;
	}

	for (i = 0; i < n; i++)
		R[i] = A[i] + B[i];
}

void matrix_scale(int m, int n, double *A, double s, double *R) 
{
    int i;
    int entries = m * n;
    
    for (i = 0; i < entries; i++)
		R[i] = A[i] * s;
}

int matrix_invert44(const double m[16], double invOut[16])
{
    double inv[16], det;
    int i;

    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return 0;

    det = 1.0 / det;

    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return 1;
}

void matrix_product332(const double A[9], const double b[3], double r[3])
{
	r[0] = A[0] * b[0] + A[3] * b[1] + A[6] * b[2];
	r[1] = A[1] * b[0] + A[4] * b[1] + A[7] * b[2];
	r[2] = A[2] * b[0] + A[5] * b[1] + A[8] * b[2];
}


double matrix_determinant3(double *A)
{
    return 
	A[0] * (A[4] * A[8] - A[5] * A[7]) -
	A[1] * (A[3] * A[8] - A[5] * A[6]) +
	A[2] * (A[3] * A[7] - A[4] * A[6]);
}



/* Create the 3x3 cross product matrix from a 3-vector */
void matrix_cross_matrix(double *v, double *v_cross) 
{
    v_cross[0] = 0.0;   v_cross[1] = -v[2]; v_cross[2] = v[1];
    v_cross[3] = v[2];  v_cross[4] = 0.0;   v_cross[5] = -v[0];
    v_cross[6] = -v[1]; v_cross[7] = v[0];  v_cross[8] = 0.0;
}

/* Cross three 4x1 vectors */
void matrix_cross4(const double *u, const double *v, const double *w, double *x) 
{
    double sub1[9] = 
	{ u[1], u[2], u[3],
	  v[1], v[2], v[3],
	  w[1], w[2], w[3] };
    	
    double sub2[9] = 
	{ u[0], u[2], u[3],
	  v[0], v[2], v[3],
	  w[0], w[2], w[3] };

    double sub3[9] = 
	{ u[0], u[1], u[3],
	  v[0], v[1], v[3],
	  w[0], w[1], w[3] };

    double sub4[9] = 
	{ u[0], u[1], u[2],
	  v[0], v[1], v[2],
	  w[0], w[1], w[2] };

    double det1 = matrix_determinant3(sub1);
    double det2 = matrix_determinant3(sub2);
    double det3 = matrix_determinant3(sub3);
    double det4 = matrix_determinant3(sub4);

    x[0] = det1;
    x[1] = det2;
    x[2] = det3;
    x[3] = det4;
}

void matrix_to_axis_angle(double *R, double *axis, double *angle)
{
    double d1 = R[7] - R[5];
    double d2 = R[2] - R[6];
    double d3 = R[3] - R[1];

    double norm = sqrt(d1 * d1 + d2 * d2 + d3 * d3);
    double x = (R[7] - R[5]) / norm;
    double y = (R[2] - R[6]) / norm;
    double z = (R[3] - R[1]) / norm;

    *angle = acos((R[0] + R[4] + R[8] - 1.0) * 0.5);

    axis[0] = x;
    axis[1] = y;
    axis[2] = z;
}

void axis_angle_to_matrix(double *axis, double angle, double *R)
{
    double ident[9];
    double n[9] = { 0.0, -axis[2], axis[1],
		    axis[2], 0.0, -axis[0],
		    -axis[1], axis[0], 0.0 };

    double nsq[9], sn[9], cnsq[9], tmp[9];

    double c, s;
    
    c = cos(angle);
    s = sin(angle);

    matrix_ident(3, ident);
    matrix_product33(n, n, nsq);
    matrix_scale(3, 3, n, s, sn);
    matrix_scale(3, 3, nsq, (1 - c), cnsq);

    matrix_sum(3, 3, 3, 3, ident, sn, tmp);
    matrix_sum(3, 3, 3, 3, tmp, cnsq, R);
}

void matrix_to_quaternion(double *R, double *q) 
{
    double n4; // the norm of quaternion multiplied by 4 
    double tr = R[0] + R[4] + R[8]; // trace of martix
    double factor;

    if (tr > 0.0) {
        q[1] = R[5] - R[7];
        q[2] = R[6] - R[2];
        q[3] = R[1] - R[3];
        q[0] = tr + 1.0;
        n4 = q[0];
    } else if ((R[0] > R[4]) && (R[0] > R[8])) {
        q[1] = 1.0 + R[0] - R[4] - R[8];
        q[2] = R[3] + R[1];
        q[3] = R[6] + R[2];
        q[0] = R[5] - R[7];
        n4 = q[1];
    } else if (R[4] > R[8]) {
        q[1] = R[3] + R[1];
        q[2] = 1.0 + R[4] - R[0] - R[8];
        q[3] = R[7] + R[5];
        q[0] = R[6] - R[2]; 
        n4 = q[2];
    } else {
        q[1] = R[6] + R[2];
        q[2] = R[7] + R[5];
        q[3] = 1.0 + R[8] - R[0] - R[4];
        q[0] = R[1] - R[3];
        n4 = q[3];
    }

	factor = 0.5 / sqrt(n4);
    q[0] *= factor;
    q[1] *= factor;
    q[2] *= factor;
    q[3] *= factor;
}

void quaternion_to_matrix(double *q, double *R)
{
    double sqw = q[0] * q[0];
    double sqx = q[1] * q[1]; 
    double sqy = q[2] * q[2];
    double sqz = q[3] * q[3];

    double tmp1, tmp2;

    R[0] =  sqx - sqy - sqz + sqw; // since sqw + sqx + sqy + sqz = 1
    R[4] = -sqx + sqy - sqz + sqw;    
    R[8] = -sqx - sqy + sqz + sqw;

    tmp1 = q[1] * q[2];
    tmp2 = q[3] * q[0];

    R[1] = 2.0 * (tmp1 + tmp2);    
    R[3] = 2.0 * (tmp1 - tmp2);        

    tmp1 = q[1] * q[3];    
    tmp2 = q[2] * q[0];    

    R[2] = 2.0 * (tmp1 - tmp2);    
    R[6] = 2.0 * (tmp1 + tmp2);    

    tmp1 = q[2] * q[3];    
    tmp2 = q[1] * q[0];    

    R[5] = 2.0 * (tmp1 + tmp2);   
    R[7] = 2.0 * (tmp1 - tmp2);
}



