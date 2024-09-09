#ifndef foster_helper_h
#define foster_helper_h
#include <iostream>
#include <cmath>
#include <fstream>
#include "json.hpp"
#include "LUDv.h"


using namespace std;
using json = nlohmann::json;

struct Atmosphere
{
    double geopotential_altitude;
    double temperature;
    double pressure;
    double density;
    double speed_of_sound;
    double dynamic_viscosity;
};

const double pi = 3.14159265358979323846264338327950288419716939937510;

// Atmosphere Functions
void get_atmospheric_properties_si(double altitude, Atmosphere& atm);
void get_atmospheric_properties_english(double altitude, Atmosphere& atm);
double gravity_si(double altitude);
double gravity_english(double altitude);
void quat_mult(double* A, double* B, double* ans);
void quat_norm(double* quat);
void euler_to_quat(double* eul, double* quat);
void quat_to_euler(double* quat, double* eul);
void body_to_fixed(double* vec, double* quat, double* ans);
void fixed_to_body(double* vec, double* quat, double* ans);
void array_copy(double* A, double* B, int size);
bool matrix_invert(double a[3][3], double b[3][3]);
void matrix_vector_mult_3(double rm[3][3], double v[3], double* ans);
void vector_cross_3(double a[3], double b[3], double ans[3]);
void print_array(double* vec, int size);
void matrix_AxB_solve(double** A, double* B, int size, double* x);
void array_copy(double* A, double* B, int size);
bool key_exists(json j_object, string key);


#endif