#include "foster_helper.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>



void get_atmospheric_properties_si(double altitude, Atmosphere& atm)
{
	double Z,T,p,rho,a,mu;

	//cosntants
	double R_E = 6356766;
	double R = 287.0528;
	double g_0 = 9.806645;
	double p_0 = 101325;

	//geopotential[m]
	Z = (R_E * altitude) / (R_E + altitude);

	//Temperature[K] and pressure
	if (Z < 11000) 
	{
		T = 288.150 + -0.0065 * (Z);
		p = p_0 * pow(((288.150 + -0.0065 * (Z)) / 288.150), -1 * g_0 / (R * -0.0065));
	} 
	else if ( Z < 20000) 
	{
		T = 216.650;
		p = p_0 * pow(((288.150 + -0.0065 * (11000)) / 288.150), -1 * g_0 / (R * -0.0065));
		p = p * exp((-1 * g_0 * (Z - 11000)) / (R * 216.650));
	} 
	else if (Z < 32000) 
	{
		T = 216.650 + 0.001 * (Z-20000);
		p = p_0 * pow(((288.150 + -0.0065 * (11000)) / 288.150), -1 * g_0 / (R * -0.0065));
		p = p * exp((-1 * g_0 * (20000 - 11000)) / (R * 216.650));
		p = p * pow(((216.65 + 0.001 * (Z - 20000)) / 216.65), -1 * g_0 / (R * 0.001));

	} 
	else if (Z < 47000) 
	{
		T = 228.650 + 0.0028 * (Z - 32000);
		p = p_0 * pow(((288.150 + -0.0065 * (11000)) / 288.150), -1 * g_0 / (R * -0.0065));
		p = p * exp((-1 * g_0 * (20000 - 11000)) / (R * 216.650));
		p = p * pow(((216.65 + 0.001 * (32000 - 20000)) / 216.65), -1 * g_0 / (R * 0.001));
		p = p * pow(((228.65 + 0.0028 * (Z - 32000)) / 228.65), -1 * g_0 / (R * 0.0028));
	}
	else if (Z < 52000) 
	{
		T = 270.65;
		p = p_0 * pow(((288.150 + -0.0065 * (11000)) / 288.150), -1 * g_0 / (R * -0.0065));
		p = p * exp((-1 * g_0 * (20000 - 11000)) / (R * 216.650));
		p = p * pow(((216.65 + 0.001 * (32000 - 20000)) / 216.65), -1 * g_0 / (R * 0.001));
		p = p * pow(((228.65 + 0.0028 * (47000 - 32000)) / 228.65), -1 * g_0 / (R * 0.0028));
		p = p * exp((-1 * g_0 * (Z - 47000)) / (R * 270.65));
	}
	else if (Z < 61000)
	{
		T = 270.65 + -0.002 * (Z - 52000);
		p = p_0 * pow(((288.150 + -0.0065 * (11000)) / 288.150), -1 * g_0 / (R * -0.0065));
		p = p * exp((-1 * g_0 * (20000 - 11000)) / (R * 216.650));
		p = p * pow(((216.65 + 0.001 * (32000 - 20000)) / 216.65), -1 * g_0 / (R * 0.001));
		p = p * pow(((228.65 + 0.0028 * (47000 - 32000)) / 228.65), -1 * g_0 / (R * 0.0028));
		p = p * exp((-1 * g_0 * (52000 - 47000)) / (R * 270.65));
		p = p * pow(((270.65 - 0.002 * (Z - 52000)) / 270.65), g_0 / (R * 0.002));
	}
	else if (Z < 79000)
	{
		T = 252.65 + -0.004 * (Z - 61000);
		p = p_0 * pow(((288.150 + -0.0065 * (11000)) / 288.150), -1 * g_0 / (R * -0.0065));
		p = p * exp((-1 * g_0 * (20000 - 11000)) / (R * 216.650));
		p = p * pow(((216.65 + 0.001 * (32000 - 20000)) / 216.65), -1 * g_0 / (R * 0.001));
		p = p * pow(((228.65 + 0.0028 * (47000 - 32000)) / 228.65), -1 * g_0 / (R * 0.0028));
		p = p * exp((-1 * g_0 * (52000 - 47000)) / (R * 270.65));
		p = p * pow(((270.65 + -0.002 * (61000 - 52000)) / 270.65), -1 * g_0 / (R * -0.002));
		p = p * pow(((252.65 + -0.004 * (Z - 61000)) / 252.65), -1 * g_0 / (R * -0.004));
	}
	else
	{
		T = 180.65;
		p = p_0 * pow(((288.150 + -0.0065 * (11000)) / 288.150), -1 * g_0 / (R * -0.0065));
		p = p * exp((-1 * g_0 * (20000 - 11000)) / (R * 216.650));
		p = p * pow(((216.65 + 0.001 * (32000 - 20000)) / 216.65), -1 * g_0 / (R * 0.001));
		p = p * pow(((228.65 + 0.0028 * (47000 - 32000)) / 228.65), -1 * g_0 / (R * 0.0028));
		p = p * exp((-1 * g_0 * (52000 - 47000)) / (R * 270.65));
		p = p * pow(((270.65 + -0.002 * (61000 - 52000)) / 270.65), -1 * g_0 / (R * -0.002));
		p = p * pow(((252.65 + -0.004 * (79000 - 61000)) / 252.65), -1 * g_0 / (R * -0.004));
		p = p * exp((-1 * g_0 * (Z - 79000)) / (R * 180.650));

	}

	// Density
	rho = p / (R * T);

	//speed of sound
	a = sqrt(1.4 * R * T);

	//dynamic viscosity
	double T0 = 273.15;
	double mu0 = 0.00001716;
	double C = 110.4;
	mu = mu0*(T0 + C)/(T + C)*pow((T/T0), 1.5);

	atm.geopotential_altitude = Z;
	atm.temperature = T;
	atm.pressure = p;
	atm.density = rho;
	atm.speed_of_sound = a;
	atm.dynamic_viscosity = mu;

}

void get_atmospheric_properties_english(double altitude, Atmosphere& atm)
{
    altitude *= 0.3048;
    get_atmospheric_properties_si(altitude, atm);
    atm.geopotential_altitude /= 0.3048;
    atm.temperature *= 1.8;
    atm.pressure *= 0.020885434304801722;
    atm.density *= 0.00194032032363104;
    atm.speed_of_sound /= 0.3048;
	atm.dynamic_viscosity /= 47.88025898;
}


double gravity_si(double H){
	double g = 9.806645 * pow((6356766 / (6356766 + H)), 2);
	return g;
}

double gravity_english(double H){
	H *= 0.3048;
	double g = 9.806645 * pow((6356766 / (6356766 + H)), 2);
	g /= 0.3048;
	return g;
}

void quat_mult(double* A, double* B, double* ans){

	ans[0] = A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3];
	ans[1] = A[0]*B[1] + A[1]*B[0] + A[2]*B[3] - A[3]*B[2]; 
	ans[2] = A[0]*B[2] - A[1]*B[3] + A[2]*B[0] + A[3]*B[1];
	ans[3] = A[0]*B[3] + A[1]*B[2] - A[2]*B[1] + A[3]*B[0];

}

void quat_norm(double* quat){
	
	double magnitude = sqrt(pow(quat[0], 2) + pow(quat[1], 2) + pow(quat[2], 2) + pow(quat[3], 2));
	quat[0] /= magnitude;
	quat[1] /= magnitude;
	quat[2] /= magnitude;
	quat[3] /= magnitude;
}

void euler_to_quat(double* eul, double* quat){

	double cphi = cos(eul[0]/2);
	double ctheta = cos(eul[1]/2);
	double cpsi = cos(eul[2]/2);
	double sphi = sin(eul[0]/2);
	double stheta = sin(eul[1]/2);
	double spsi =  sin(eul[2]/2);
	
	quat[0] = cphi * ctheta * cpsi + sphi * stheta * spsi;
	quat[1] = sphi * ctheta * cpsi - cphi * stheta * spsi;
	quat[2] = cphi * stheta * cpsi + sphi * ctheta * spsi;
	quat[3] = cphi * ctheta * spsi - sphi * stheta * cpsi;
}

void quat_to_euler(double* quat, double* eul){


	if (quat[0]*quat[2] - quat[1]*quat[3] == 0.5){
		eul[0] = 2 * asin(quat[1] / cos(pi / 4.));
		eul[1] = pi / 2.0;
		eul[2] = 0.0;
	}
	else if (quat[0]*quat[2] - quat[1]*quat[3] == -0.5)
	{
		eul[0] = -2*asin(quat[1]/cos(pi / 4.));
		eul[1] = -pi/2.0;
		eul[2] = 0.0;
	}
	else
	{
		eul[0] = atan2(2 * (quat[0]*quat[1] + quat[2]*quat[3]), (pow(quat[0],2) + pow(quat[1],2) - pow(quat[2],2) - pow(quat[3],2)));
		eul[1] = asin(2 * (quat[0]*quat[2] - quat[1]*quat[3]));
		eul[2] = atan2(2 * (quat[0]*quat[3] + quat[2]*quat[1]), (pow(quat[0],2) + pow(quat[1],2) - pow(quat[2],2) - pow(quat[3], 2)));
	}
	
}

void body_to_fixed(double* vec, double* quat, double* ans){

	ans[0] = (pow(quat[1], 2) + pow(quat[0], 2) - pow(quat[2], 2) - pow(quat[3], 2)) * vec[0] + 2 * (quat[1] * quat[2] - quat[3] * quat[0]) * vec[1] + 2 * (quat[1] * quat[3] + quat[2] * quat[0]) * vec[2];
	ans[1] = 2 * (quat[1] * quat[2] + quat[3] * quat[0]) * vec[0] + (pow(quat[2], 2) + pow(quat[0], 2) - pow(quat[1], 2) - pow(quat[3], 2)) * vec[1] + 2 * (quat[2] * quat[3] - quat[1] * quat[0]) * vec[2];
	ans[2] = 2 * (quat[1] * quat[3] - quat[2] * quat[0]) * vec[0] + 2 * (quat[2] * quat[3] + quat[1] * quat[0]) * vec[1] + (pow(quat[3], 2) + pow(quat[0], 2) - pow(quat[1], 2) - pow(quat[2], 2)) * vec[2];
}

void fixed_to_body(double* vec, double* quat, double* ans){
	
	ans[0] = (pow(quat[1], 2) + pow(quat[0], 2) - pow(quat[2], 2) - pow(quat[3], 2)) * vec[0] + 2 * (quat[1] * quat[2] + quat[3] * quat[0]) * vec[1] + 2 * (quat[1] * quat[3] - quat[2] * quat[0]) * vec[2];
	ans[1] = 2 * (quat[1] * quat[2] - quat[3] * quat[0]) * vec[0] + (pow(quat[2], 2) + pow(quat[0], 2) - pow(quat[1], 2) - pow(quat[3], 2)) * vec[1] + 2 * (quat[2] * quat[3] + quat[1] * quat[0]) * vec[2];
	ans[2] = 2 * (quat[1] * quat[3] + quat[2] * quat[0]) * vec[0] + 2 * (quat[2] * quat[3] - quat[1] * quat[0]) * vec[1] + (pow(quat[3], 2) + pow(quat[0], 2) - pow(quat[1], 2) - pow(quat[2], 2)) * vec[2];
}

bool matrix_invert(double a[3][3], double b[3][3]) {
	double det = a[0][0]*(a[1][1]*a[2][2] - a[1][2]*a[2][1])
	- a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0])
	+ a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0]);

	if (det == 0) {
		// The matrix is singular, cannot be inverted
		return false;
	}

	b[0][0] = (a[1][1]*a[2][2] - a[1][2]*a[2][1]) / det;
	b[0][1] = (a[0][2]*a[2][1] - a[0][1]*a[2][2]) / det;
	b[0][2] = (a[0][1]*a[1][2] - a[0][2]*a[1][1]) / det;
	b[1][0] = (a[1][2]*a[2][0] - a[1][0]*a[2][2]) / det;
	b[1][1] = (a[0][0]*a[2][2] - a[0][2]*a[2][0]) / det;
	b[1][2] = (a[0][0]*a[1][2] - a[0][2]*a[1][0]) / det;
	b[2][0] = (a[1][0]*a[2][1] - a[1][1]*a[2][0]) / det;
	b[2][1] = (a[0][1]*a[2][0] - a[0][0]*a[2][1]) / det;
	b[2][2] = (a[0][0]*a[1][1] - a[0][1]*a[1][0]) / det;

return true;
}

void matrix_vector_mult_3(double rm[3][3], double v[3], double* ans){
	int i, j;
	for (i = 0; i < 3; i++){
		ans[i] = 0;
		for (j = 0; j < 3; j++){
			ans[i] = ans[i] + rm[i][j]*v[j];
		}
	}
}

void vector_cross_3(double a[3], double b[3], double ans[3]){

ans[0] = a[1]*b[2] - a[2]*b[1];
ans[1] = a[2]*b[0] - a[0]*b[2];
ans[2] = a[0]*b[1] - a[1]*b[0];
}

void print_array(double* vec, int size){

	int i;
	for (i=0; i< size; i++){
		printf("%20.12e, " , vec[i]);
		
	}
	cout<<endl;
}

void matrix_AxB_solve(double** A, double* B, int size, double* x){
	double* Pvt = new double[size];
	double* temp = new double[size];
	array_copy(B,temp,size);
	LUDecomp(A, size, Pvt);
	LUSolve( A, size, B, Pvt);
	array_copy(B,x,size);
	array_copy(temp,B,size);
}

void array_copy(double* A, double* B, int size){

	for (int i = 0; i < size; i++) {B[i] = A[i];}
}

bool key_exists(json j_object, string key){
return (j_object.find(key) != j_object.end());
}