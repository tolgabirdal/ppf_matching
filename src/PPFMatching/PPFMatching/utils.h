
#ifndef __UTILS_H_
#define __UTILS_H_

// Useful Macros
#define TNORM3(v) (sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]))

#define TNormalize3(v)\
double normTemp=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);\
if (normTemp>0)\
{\
v[0]/=normTemp;\
v[1]/=normTemp;\
v[2]/=normTemp;\
}

#if defined (__cplusplus)
extern "C" {
#endif 

	void matrix_cross4(const double *u, const double *v, const double *w, double *x);
	void matrix_cross_matrix(double *v, double *v_cross) ;
	
	void matrix_to_axis_angle(double *R, double *axis, double *angle);
	void axis_angle_to_matrix(double *axis, double angle, double *R);
	void matrix_to_quaternion(double *R, double *q);
	void quaternion_to_matrix(double *q, double *R);

	static __inline void matrix_product33(double *A, double *B, double *R)
	{
		R[0] = A[0] * B[0] + A[1] * B[3] + A[2] * B[6];
		R[1] = A[0] * B[1] + A[1] * B[4] + A[2] * B[7];
		R[2] = A[0] * B[2] + A[1] * B[5] + A[2] * B[8];    

		R[3] = A[3] * B[0] + A[4] * B[3] + A[5] * B[6];
		R[4] = A[3] * B[1] + A[4] * B[4] + A[5] * B[7];
		R[5] = A[3] * B[2] + A[4] * B[5] + A[5] * B[8];    

		R[6] = A[6] * B[0] + A[7] * B[3] + A[8] * B[6];
		R[7] = A[6] * B[1] + A[7] * B[4] + A[8] * B[7];
		R[8] = A[6] * B[2] + A[7] * B[5] + A[8] * B[8];            
	}

#if defined (__cplusplus)
}
#endif 

#endif