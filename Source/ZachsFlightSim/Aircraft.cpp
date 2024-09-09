// Fill out your copyright notice in the Description page of Project Settings.


#include "Aircraft.h"

// Sets default values
AAircraft::AAircraft()
{
 	// Set this pawn to call Tick() every frame.  You can turn this off to improve performance if you don't need it.
	PrimaryActorTick.bCanEverTick = true;

}

// Called when the game starts or when spawned
void AAircraft::BeginPlay()
{
	Super::BeginPlay();
	
}

// Called every frame
void AAircraft::Tick(float DeltaTime)
{
	Super::Tick(DeltaTime);

}

// Called to bind functionality to input
void AAircraft::SetupPlayerInputComponent(UInputComponent* PlayerInputComponent)
{
	Super::SetupPlayerInputComponent(PlayerInputComponent);

}

void AAircraft::rk4Solver(double t, double* y0, double dt, int size, double* FM, double* ans)
{

    // calculate k1
    PhysicsModel(t, y0, FM, m_k1);
    //cout<<"mk1"<<endl;
    //print_array(m_k1, m_size);
    //std::cout <<scientific << m_k1[2] << std::endl;

    for (int i = 0; i < size; i++) {
        m_k1[i] = dt * m_k1[i];
    }

    // calculate y + k1/2
    for (int i = 0; i < size; i++) {
        ans[i] = y0[i] + m_k1[i] / 2;
    }

    // calculate k2
    PhysicsModel(t, ans, FM, m_k2);
    //cout<<"mk2"<<endl;
    //print_array(m_k2, m_size);
    //std::cout <<scientific<< m_k2[2] << std::endl;

    for (int i = 0; i < size; i++) {
        m_k2[i] = dt * m_k2[i];
    }

    // calculate y + k2/2
    for (int i = 0; i < size; i++) {
        ans[i] = y0[i] + m_k2[i] / 2;
    }

    // calculate k3
    PhysicsModel(t, ans, FM, m_k3);
    //cout<<"mk3"<<endl;
    //print_array(m_k3, m_size);
    //std::cout <<scientific<< m_k3[2] << std::endl;

    for (int i = 0; i < size; i++) {
        m_k3[i] = dt * m_k3[i];
    }

    // calculate y + k3
    for (int i = 0; i < size; i++) {
        ans[i] = y0[i] + m_k3[i];
    }

    // calculate k4
    PhysicsModel(t, ans, FM, m_k4);
    //cout<<"mk4"<<endl;
    //print_array(m_k4, m_size);
    //std::cout <<scientific<< m_k4[2] << std::endl;

    for (int i = 0; i < size; i++) {
        m_k4[i] = dt * m_k4[i];
    }

    // calculate ans
    for (int i = 0; i < size; i++) {
        ans[i] = y0[i] + (m_k1[i] / 6 + m_k2[i] / 3 + m_k3[i] / 3 + m_k4[i] / 6);
    }

}
void AAircraft::PhysicsModel(double t, double* y0, double* FM, double* ans)
{

    aerodynamicModel(y0, FM);

    double u = y0[0];
    double v = y0[1];
    double w = y0[2];
    double p = y0[3];
    double q = y0[4];
    double r = y0[5];
    double xf = y0[6];
    double yf = y0[7];
    double zf = y0[8];
    double e0 = y0[9];
    double ex = y0[10];
    double ey = y0[11];
    double ez = y0[12];

    ans[0] = gravity_english(-zf) / m_weight * (FM[0]) + 2 * (ex * ez - ey * e0) * gravity_english(-zf) + (r * v - q * w);
    ans[1] = gravity_english(-zf) / m_weight * (FM[1]) + 2 * (ey * ez + ex * e0) * gravity_english(-zf) + (p * w - r * u);
    ans[2] = gravity_english(-zf) / m_weight * (FM[2]) + (ez * ez + e0 * e0 - ex * ex - ey * ey) * gravity_english(-zf) + (q * u - p * v);


    double pqrTemp[3];
    pqrTemp[0] = FM[3] + (m_Iyy - m_Izz) * q * r + m_Iyz * (q * q - r * r) + m_Ixz * p * q - m_Ixy * p * r;
    pqrTemp[1] = FM[4] + (m_Izz - m_Ixx) * p * r + m_Ixz * (r * r - p * p) + m_Ixy * q * r - m_Iyz * p * q;
    pqrTemp[2] = FM[5] + (m_Ixx - m_Iyy) * p * q + m_Ixy * (p * p - q * q) + m_Iyz * p * r - m_Ixz * q * r;

    double hpqrTemp[3];
    matrix_vector_mult_3(m_hmatrix, &y0[3], hpqrTemp);
    hpqrTemp[0] += pqrTemp[0];
    hpqrTemp[1] += pqrTemp[1];
    hpqrTemp[2] += pqrTemp[2];

    matrix_vector_mult_3(m_ImatrixInv, hpqrTemp, &ans[3]);

    double A[4];
    A[0] = 0.0;
    A[1] = u;
    A[2] = v;
    A[3] = w;
    double B[4];
    B[0] = e0;
    B[1] = -ex;
    B[2] = -ey;
    B[3] = -ez;
    double C[4];
    double D[4];
    quat_mult(A, B, C);
    quat_mult(&y0[9], C, D);

    ans[6] = D[1];
    ans[7] = D[2];
    ans[8] = D[3];
    ans[9] = 0.5 * (-ex * p - ey * q - ez * r);
    ans[10] = 0.5 * (e0 * p - ez * q + ey * r);
    ans[11] = 0.5 * (ez * p + e0 * q - ex * r);
    ans[12] = 0.5 * (-ey * p + ex * q + e0 * r);

    //UE_LOG(LogTemp, Warning, TEXT("forces and moments = %f\t%f\t%f"), ans[3], ans[4], ans[5]);

}
void AAircraft::aerodynamicModel(double* y0, double* FM)
{

    get_atmospheric_properties_english(-y[8], m_atm);


    //angle of atack and velocity
    double V = sqrt(y0[0] * y0[0] + y0[1] * y0[1] + y0[2] * y0[2]);
    double alpha = atan2(y0[2], y0[0]);
    double beta = asin(y0[1] / V);

    //nondimentionalize p q and r
    double pbar = (y0[3] * 0.5 * m_wing_span) / V;
    double qbar = (y0[4] * 0.5 * m_wing_chord) / V;
    double rbar = (y0[5] * 0.5 * m_wing_span) / V;

    //calculate aerodynamic coefficents below stall
    double CL1 = m_CL0 + m_CLa * alpha;
    double CL = m_CL0 + m_CLa * alpha + m_CLqbar * qbar + m_CLde * m_control[0];
    double CS = m_CSb * beta + m_CSpbar * pbar + m_CSrbar * rbar + m_CSda * m_control[1] + m_CSdr * m_control[2];
    double CD = m_CDL0 + m_CDL * CL1 + m_CDL2 * CL1 * CL1 + m_CDS2 * CS * CS + (m_CDLqbar * CL1 + m_CDqbar) * qbar + (m_CDLde * CL1 + m_CDde) * m_control[0] + m_CDde2 * m_control[0] * m_control[0];
    double Cl = m_Clb * beta + m_Clpbar * pbar + (m_ClLrbar * CL1 + m_Clrbar) * rbar + m_Clda * m_control[1] + m_Cldr * m_control[2];
    double Cm = m_Cm0 + m_Cma * alpha + m_Cmqbar * qbar + m_Cmde * m_control[0];
    double Cn = m_Cnb * beta + (m_CnLpbar * CL1 + m_Cnpbar) * pbar + m_Cnrbar * rbar + (m_CnLda * CL1 + m_Cnda) * m_control[1] + m_Cndr * m_control[2];
    
    //calculate aerodynamic coefficents for a flat plate
    double sign;

    if (alpha>0) {
        sign = 1;
    }
    else {
        sign = -1;
    }

    double CLfp = 2 * sign * sin(alpha) * sin(alpha) * cos(alpha);
    double CDfp = 2 * pow(sin(abs(alpha)), 3 / 2);
    double Cmfp = -0.8 * sin(alpha);

    //calculate CL, CD, Cm including stall model.

    //CL
    double M;
    double ad;

    M = 6.807527802;
    ad = 0.691289282;


    double sigma = (1 + exp(-M * (alpha - ad)) + exp(M * (alpha + ad))) / ((1 + exp(-M * (alpha - ad))) * (1 + exp(M * (alpha + ad))));
    
    CL = (1 - sigma) * CL + sigma * CLfp;

    //CD
    M = 2.301926314;
    ad = 0.038013541;

    sigma = (1 + exp(-M * (alpha - ad)) + exp(M * (alpha + ad))) / ((1 + exp(-M * (alpha - ad))) * (1 + exp(M * (alpha + ad))));

    CD = (1 - sigma) * CD + sigma * CDfp;

    //Cm
    M = 4.197575289;
    ad = 1.275857798;

    sigma = (1 + exp(-M * (alpha - ad)) + exp(M * (alpha + ad))) / ((1 + exp(-M * (alpha - ad))) * (1 + exp(M * (alpha + ad))));

    Cm = (1 - sigma) * Cm + sigma * Cmfp;

    // Compressablity corrections

    double sonicVel = m_atm.speed_of_sound;
    double Mach = V / sonicVel;

    CL = CL / (pow(1 - Mach * Mach, 0.5) + ((Mach * Mach) * (1 + .4 / 2 * Mach * Mach)) / (1 + pow(1 - Mach * Mach, 0.5)) * CL / 2);
    CD = CD / (pow(1 - Mach * Mach, 0.5) + ((Mach * Mach) * (1 + .4 / 2 * Mach * Mach)) / (1 + pow(1 - Mach * Mach, 0.5)) * CD / 2);
    Cm = Cm / (pow(1 - Mach * Mach, 0.5) + ((Mach * Mach) * (1 + .4 / 2 * Mach * Mach)) / (1 + pow(1 - Mach * Mach, 0.5)) * Cm / 2);

    // Calculate thrust forces
    double thrust[6];
    double thrust_location[3];
    thrust_location[0] = m_thrust_location_x;
    thrust_location[1] = m_thrust_location_y;
    thrust_location[2] = m_thrust_location_z;
    double total_thrust = m_control[3] * pow((m_atm.density / m_rho0), m_control[4]) * m_control[5];
    thrust[0] = m_thrust_direction_x * total_thrust;
    thrust[1] = m_thrust_direction_y * total_thrust;
    thrust[2] = m_thrust_direction_z * total_thrust;
    //calculates thrust moments
    vector_cross_3(thrust_location, thrust, &thrust[3]);

    //CG Shift
    double CG_shift[3];
    CG_shift[0] = m_CG_shift_x;
    CG_shift[1] = m_CG_shift_y;
    CG_shift[2] = m_CG_shift_z;

    double constant = 0.5 * m_atm.density * V * V * m_wing_area;

    FM[0] = thrust[0] + constant * (CL * sin(alpha) - CS * cos(alpha) * sin(beta) - CD * cos(alpha) * cos(beta));
    FM[1] = thrust[1] + constant * (CS * cos(beta) - CD * sin(beta));
    FM[2] = thrust[2] + constant * (-CL * cos(alpha) - CS * sin(alpha) * sin(beta) - CD * sin(alpha) * cos(beta));

    //CG shift
    double CG_ShiftCross[3];
    vector_cross_3(CG_shift, FM, CG_ShiftCross);

    FM[3] = (thrust[3] + constant * (m_wing_span * Cl)) - CG_ShiftCross[0];
    FM[4] = (thrust[4] + constant * (m_wing_chord * Cm)) - CG_ShiftCross[1];
    FM[5] = (thrust[5] + constant * (m_wing_span * Cn)) - CG_ShiftCross[2];

    //UE_LOG(LogTemp, Warning, TEXT("moment coef = %f\t%f\t%f"), FM[3], FM[4], FM[5]);
}
void AAircraft::initalizeState()
{

    double elevation_angle = data["initial"]["state"]["elevation_angle[deg]"];
    double bank_angle = data["initial"]["state"]["bank_angle[deg]"];
    double alpha = data["initial"]["state"]["alpha[deg]"];
    double beta = data["initial"]["state"]["beta[deg]"];
    double p = data["initial"]["state"]["p[deg/s]"];
    double q = data["initial"]["state"]["q[deg/s]"];
    double r = data["initial"]["state"]["r[deg/s]"];
    double ailoron = data["initial"]["state"]["aileron[deg]"];
    double elevator = data["initial"]["state"]["elevator[deg]"];
    double rudder = data["initial"]["state"]["rudder[deg]"];
    double throttle = data["initial"]["state"]["throttle"];

    elevation_angle *= pi / 180;
    bank_angle *= pi / 180;
    m_heading_angle *= pi / 180;
    alpha *= pi / 180;
    beta *= pi / 180;
    p *= pi / 180;
    q *= pi / 180;
    r *= pi / 180;
    ailoron *= pi / 180;
    elevator *= pi / 180;
    rudder *= pi / 180;

    // initial state
    y[0] = m_airspeed * cos(alpha) * cos(beta);
    y[1] = m_airspeed * sin(beta);
    y[2] = m_airspeed * sin(alpha) * cos(beta);
    y[3] = p;
    y[4] = q;
    y[5] = r;
    y[6] = 0.0;
    y[7] = 0.0;
    y[8] = -1 * m_altitude;
    y[9] = bank_angle;
    y[10] = elevation_angle;
    y[11] = m_heading_angle;
    y[12] = 0.0;

    //cout<<"init"<<endl;
    //cout<<init[0]<<endl;

    euler_to_quat(&y[9], &y[9]);

    //m_control vector
    m_control[0] = elevator;
    m_control[1] = ailoron;
    m_control[2] = rudder;
    m_control[3] = throttle;
    m_control[4] = m_a;
    m_control[5] = m_T0;

    //inertial tensor
    m_Imatrix[0][0] = m_Ixx;
    m_Imatrix[0][1] = -m_Ixy;
    m_Imatrix[0][2] = -m_Ixz;
    m_Imatrix[1][0] = -m_Ixy;
    m_Imatrix[1][1] = m_Iyy;
    m_Imatrix[1][2] = -m_Iyz;
    m_Imatrix[2][0] = -m_Ixz;
    m_Imatrix[2][1] = -m_Iyz;
    m_Imatrix[2][2] = m_Izz;

    matrix_invert(m_Imatrix, m_ImatrixInv);

    //gyrocopic tensor
    m_hmatrix[0][0] = 0.0;
    m_hmatrix[0][1] = -m_hz;
    m_hmatrix[0][2] = m_hy;
    m_hmatrix[1][0] = m_hz;
    m_hmatrix[1][1] = 0.0;
    m_hmatrix[1][2] = -m_hx;
    m_hmatrix[2][0] = -m_hy;
    m_hmatrix[2][1] = m_hx;
    m_hmatrix[2][2] = 0.0;
}
void AAircraft::initalizeTrim()
{
    //inertial tensor
    m_Imatrix[0][0] = m_Ixx;
    m_Imatrix[0][1] = -m_Ixy;
    m_Imatrix[0][2] = -m_Ixz;
    m_Imatrix[1][0] = -m_Ixy;
    m_Imatrix[1][1] = m_Iyy;
    m_Imatrix[1][2] = -m_Iyz;
    m_Imatrix[2][0] = -m_Ixz;
    m_Imatrix[2][1] = -m_Iyz;
    m_Imatrix[2][2] = m_Izz;

    matrix_invert(m_Imatrix, m_ImatrixInv);

    //gyrocopic tensor
    m_hmatrix[0][0] = 0.0;
    m_hmatrix[0][1] = -m_hz;
    m_hmatrix[0][2] = m_hy;
    m_hmatrix[1][0] = m_hz;
    m_hmatrix[1][1] = 0.0;
    m_hmatrix[1][2] = -m_hx;
    m_hmatrix[2][0] = -m_hy;
    m_hmatrix[2][1] = m_hx;
    m_hmatrix[2][2] = 0.0;


    // Solver parameters
    double delta_newton = data["initial"]["trim"]["solver"]["finite_difference_step_size"];
    double relaxation_factor = data["initial"]["trim"]["solver"]["relaxation_factor"];
    double tolerance = data["initial"]["trim"]["solver"]["tolerance"];
    int max_iterations = data["initial"]["trim"]["solver"]["max_iterations"];
    bool verbose = data["initial"]["trim"]["solver"]["verbose"];

    m_heading_angle *= pi / 180;

    double alpha = 0;
    double beta = 0;
    m_control[0] = 0;
    m_control[1] = 0;
    m_control[2] = 0;
    m_control[3] = 0;
    m_control[4] = m_a;
    m_control[5] = m_T0;

    double p = 0;
    double q = 0;
    double r = 0;
    double u = 0;
    double v = 0;
    double w = 0;
    double bank_angle = 0;
    double RotConst = 0;
    int iter = 0;
    double elevation_angle1 = 0;
    double elevation_angle2 = 0;
    double check1 = 0;
    double check2 = 0;


    // Calculates constants and creates variables
    double g = gravity_english(m_altitude);
    get_atmospheric_properties_english(m_altitude, m_atm);

    double rho = m_atm.density;
    double X, Y;
    X = 0;
    Y = 0;
    double y0[9];
    double G[6];
    double FM[6];
    double R[6];
    double Rp[6];
    double Rn[6];
    double deltaG[6];

    double** J = new double* [6];
    for (int i = 0; i < 6; i++) {
        J[i] = new double[6];
    }


    bool notConverged = true;

    double elevation_angle;
    double climb_angle;
    bool side_slip = false;



    if (key_exists(data["initial"]["trim"], "bank_angle[deg]")) {
        bank_angle = data["initial"]["trim"]["bank_angle[deg]"];
        bank_angle *= pi / 180;

    }
    else {
        beta = data["initial"]["trim"]["sideslip_angle[deg]"];
        beta *= pi / 180;
        side_slip = true;
    }


    do {
        iter++;

        // calculates body fixed airspeed
        u = m_airspeed * cos(alpha) * cos(beta);
        v = m_airspeed * sin(beta);
        w = m_airspeed * sin(alpha) * cos(beta);

        if (data["initial"]["trim"]["type"] == "shss") {
            p = 0;
            q = 0;
            r = 0;
        }
        else {
            // calculates p q r for steady coordinated turn
            RotConst = (g * sin(bank_angle) * cos(elevation_angle)) / (u * cos(elevation_angle) * cos(bank_angle) + w * sin(elevation_angle));
            p = RotConst * -sin(elevation_angle);
            q = RotConst * sin(bank_angle) * cos(elevation_angle);
            r = RotConst * cos(bank_angle) * cos(elevation_angle);
        }

        if (key_exists(data["initial"]["trim"], "elevation_angle[deg]")) {
            elevation_angle = data["initial"]["trim"]["elevation_angle[deg]"];
            elevation_angle *= pi / 180;

        }
        else {
            climb_angle = data["initial"]["trim"]["climb_angle[deg]"];
            climb_angle *= pi / 180;
            elevation_angle1 = asin((u * m_airspeed * sin(climb_angle) + (v * sin(bank_angle) + w * cos(bank_angle)) * sqrt(u * u + ((v * sin(bank_angle) + w * cos(bank_angle)) * (v * sin(bank_angle) + w * cos(bank_angle))) - (m_airspeed * m_airspeed * sin(climb_angle) * sin(climb_angle)))) / (u * u + ((v * sin(bank_angle) + w * cos(bank_angle)) * (v * sin(bank_angle) + w * cos(bank_angle)))));
            elevation_angle2 = asin((u * m_airspeed * sin(climb_angle) - (v * sin(bank_angle) + w * cos(bank_angle)) * sqrt(u * u + ((v * sin(bank_angle) + w * cos(bank_angle)) * (v * sin(bank_angle) + w * cos(bank_angle))) - (m_airspeed * m_airspeed * sin(climb_angle) * sin(climb_angle)))) / (u * u + ((v * sin(bank_angle) + w * cos(bank_angle)) * (v * sin(bank_angle) + w * cos(bank_angle)))));

            check1 = abs(u * sin(elevation_angle1) - (v * sin(bank_angle) + w * cos(bank_angle)) * cos(elevation_angle1) - m_airspeed * sin(climb_angle));
            check2 = abs(u * sin(elevation_angle2) - (v * sin(bank_angle) + w * cos(bank_angle)) * cos(elevation_angle2) - m_airspeed * sin(climb_angle));
            if (check1 <= check2) {
                elevation_angle = elevation_angle1;
            }
            else {
                elevation_angle = elevation_angle2;
            }
        }

        // Define state vector
        y0[0] = u;
        y0[1] = v;
        y0[2] = w;
        y0[3] = p;
        y0[4] = q;
        y0[5] = r;
        y0[6] = X;
        y0[7] = Y;
        y0[8] = -m_altitude;

        // create G vector
        G[0] = alpha;
        G[3] = m_control[0];
        G[2] = m_control[1];
        G[4] = m_control[2];
        G[5] = m_control[3];

        if (side_slip) {
            G[1] = bank_angle;
        }
        else {
            G[1] = beta;
        }

        if (verbose) {
            cout << "G vector" << endl;
            print_array(G, 6);

            cout << "state" << endl;
            print_array(y0, 9);
        }

        aerodynamicModel(y0, FM);

        if (verbose) {

            cout << "Forces and Moments" << endl;
            print_array(FM, 6);
        }


        // calculate R
        GetR(G, y0, g, elevation_angle, bank_angle, beta, side_slip, FM, R);


        for (int i = 0; i < 6; i++) {

            G[i] += delta_newton;
            m_control[0] = G[3];
            m_control[1] = G[2];
            m_control[2] = G[4];
            m_control[3] = G[5];
            GetR(G, y0, g, elevation_angle, bank_angle, beta, side_slip, FM, Rp);

            G[i] -= 2 * delta_newton;

            if (verbose) {
                cout << "R" << endl;

                print_array(R, 6);

                cout << "G vector" << endl;
                print_array(G, 6);

                cout << "State" << endl;
                print_array(y, 9);

                cout << "Forces and Moments" << endl;
                print_array(FM, 6);

                cout << "R postitive" << endl;
                print_array(Rp, 6);
            }

            GetR(G, y0, g, elevation_angle, bank_angle, beta, side_slip, FM, Rn);

            if (verbose) {
                cout << "G vector" << endl;
                print_array(G, 6);

                cout << "State" << endl;
                print_array(y, 9);

                cout << "Forces and Moments" << endl;
                print_array(FM, 6);

                cout << "R negitive" << endl;
                print_array(Rn, 6);

            }

            G[i] += delta_newton;

            J[0][i] = (Rp[0] - Rn[0]) / (2 * delta_newton);
            J[1][i] = (Rp[1] - Rn[1]) / (2 * delta_newton);
            J[2][i] = (Rp[2] - Rn[2]) / (2 * delta_newton);
            J[3][i] = (Rp[3] - Rn[3]) / (2 * delta_newton);
            J[4][i] = (Rp[4] - Rn[4]) / (2 * delta_newton);
            J[5][i] = (Rp[5] - Rn[5]) / (2 * delta_newton);
        }

        matrix_AxB_solve(J, R, 6, deltaG);

        if (verbose) {

            cout << "J matrix" << endl;
            for (int j = 0; j < 6; j++) {
                print_array(&J[j][0], 6);
            }
            cout << "Delta G" << endl;
            print_array(deltaG, 6);
        }

        for (int i = 0; i < 6; i++) {
            G[i] = G[i] - relaxation_factor * deltaG[i];
        }

        if (verbose) {
            cout << "new G" << endl;
            print_array(G, 6);

        }

        alpha = G[0];

        if (side_slip) {
            bank_angle = G[1];
        }
        else {
            beta = G[1];
        }
        m_control[0] = G[3];
        m_control[1] = G[2];
        m_control[2] = G[4];
        m_control[3] = G[5];

        int n = sizeof(R) / sizeof(R[0]);

        double check = 0; // Initialize the maximum absolute value

        for (int i = 0; i < n; i++) {
            double absoluteValue = abs(R[i]); // Calculate the absolute value of the current element
            if (absoluteValue > check) {
                check = absoluteValue; // Update the maximum absolute value if needed
            }
        }

        if (tolerance > abs(check)) {
            notConverged = false;
        }

        if (max_iterations < iter) {
            notConverged = false;
        }

    } while (notConverged);

    m_control[0] = G[3];
    m_control[1] = G[2];
    m_control[2] = G[4];
    m_control[3] = G[5];
    double euler[3];

    alpha = G[0];
    if (side_slip) {
        euler[0] = G[1];
        beta = beta;
    }
    else {
        beta = G[1];
        euler[9] = bank_angle;
    }
    euler[1] = elevation_angle;
    euler[2] = 0;

    double V = m_airspeed;
    y[0] = m_airspeed * cos(G[0]) * cos(beta);
    y[1] = m_airspeed * sin(beta);
    y[2] = m_airspeed * sin(G[0]) * cos(beta);
    if (data["initial"]["trim"]["type"] == "shss") {
        y[3] = 0;
        y[4] = 0;
        y[5] = 0;
    }
    else {
        RotConst = (g * sin(bank_angle) * cos(elevation_angle)) / (y[0] * cos(elevation_angle) * cos(bank_angle) + y[2] * sin(elevation_angle));
        y[3] = RotConst * (-sin(elevation_angle));
        y[4] = RotConst * sin(bank_angle) * cos(elevation_angle);
        y[5] = RotConst * cos(bank_angle) * cos(elevation_angle);
    }
    y[6] = 0;
    y[7] = 0;
    y[8] = -m_altitude;

    euler_to_quat(euler, &y[9]);
    delete(J);

}
void AAircraft::GetR(double* G, double* y0, double g, double elevation_angle, double bank_angle, double beta, bool side_slip, double* FM, double* R)
{
    // cout<<"G vector"<<endl;
    //     print_array(G,6);
    // get forces and moments

    if (side_slip) {
        bank_angle = G[1];
    }
    else {
        beta = G[1];
    }

    // Update state vecotr
    double V = m_airspeed;
    y[0] = m_airspeed * cos(G[0]) * cos(beta);
    y[1] = m_airspeed * sin(beta);
    y[2] = m_airspeed * sin(G[0]) * cos(beta);
    if (data["initial"]["trim"]["type"] == "shss") {
        y[3] = 0;
        y[4] = 0;
        y[5] = 0;
    }
    else {
        double RotConst = (g * sin(bank_angle) * cos(elevation_angle)) / (y[0] * cos(elevation_angle) * cos(bank_angle) + y[2] * sin(elevation_angle));
        y[3] = RotConst * (-sin(elevation_angle));
        y[4] = RotConst * sin(bank_angle) * cos(elevation_angle);
        y[5] = RotConst * cos(bank_angle) * cos(elevation_angle);
    }

    // Update m_control vector
    m_control[0] = G[3];
    m_control[1] = G[2];
    m_control[2] = G[4];
    m_control[3] = G[5];

    //Update Forces and moments
    aerodynamicModel(y, FM);

    // Calculate R vector
    R[0] = FM[0] - m_weight * sin(elevation_angle) + (y[5] * y[1] - y[4] * y[2]) * m_weight / g;
    R[1] = FM[1] + m_weight * sin(bank_angle) * cos(elevation_angle) + (y[3] * y[2] - y[5] * y[0]) * m_weight / g;
    R[2] = FM[2] + m_weight * cos(bank_angle) * cos(elevation_angle) + (y[4] * y[0] - y[3] * y[1]) * m_weight / g;
    R[3] = FM[3] - m_hz * y[4] + m_hy * y[5] + (m_Iyy - m_Izz) * y[4] * y[5] + m_Iyz * (y[4] * y[4] - y[5] * y[5]) + m_Ixz * y[3] * y[4] - m_Ixy * y[3] * y[5];
    R[4] = FM[4] + m_hz * y[3] - m_hx * y[5] + (m_Izz - m_Ixx) * y[3] * y[5] + m_Ixz * (y[5] * y[5] - y[3] * y[3]) + m_Ixy * y[4] * y[5] - m_Iyz * y[3] * y[4];
    R[5] = FM[5] - m_hy * y[3] + m_hx * y[4] + (m_Ixx - m_Iyy) * y[3] * y[4] + m_Ixy * (y[3] * y[3] - y[4] * y[4]) + m_Iyz * y[3] * y[5] - m_Ixz * y[4] * y[5];
}


void AAircraft::InitalizeAircraftFromJSON(FString configFileName)
{
    string filename = TCHAR_TO_UTF8(*configFileName);
    ifstream f(filename);
    data = json::parse(f);

    //Simulation
    m_time_step = data["simulation"]["time_step[s]"];
    m_total_time = data["simulation"]["total_time[s]"];

    //Aircraft
    json m_CG_shift = data["aircraft"]["CG_shift[ft]"];
    m_CG_shift_x = m_CG_shift[0];
    m_CG_shift_y = m_CG_shift[1];
    m_CG_shift_z = m_CG_shift[2];
    m_wing_area = data["aircraft"]["wing_area[ft^2]"];
    m_wing_span = data["aircraft"]["wing_span[ft]"];
    m_wing_chord = m_wing_area / m_wing_span;

    m_weight = data["aircraft"]["weight[lbf]"];
    m_Ixx = data["aircraft"]["Ixx[slug-ft^2]"];
    m_Iyy = data["aircraft"]["Iyy[slug-ft^2]"];
    m_Izz = data["aircraft"]["Izz[slug-ft^2]"];
    m_Ixy = data["aircraft"]["Ixy[slug-ft^2]"];
    m_Ixz = data["aircraft"]["Ixz[slug-ft^2]"];
    m_Iyz = data["aircraft"]["Iyz[slug-ft^2]"];
    m_hx = data["aircraft"]["hx[slug-ft^2/s]"];
    m_hy = data["aircraft"]["hy[slug-ft^2/s]"];
    m_hz = data["aircraft"]["hz[slug-ft^2/s]"];

    //Thrust
    json m_thrust_direction = data["aircraft"]["thrust"]["direction"];
    m_thrust_direction_x = m_thrust_direction[0];
    m_thrust_direction_y = m_thrust_direction[1];
    m_thrust_direction_z = m_thrust_direction[2];
    json m_thrust_location = data["aircraft"]["thrust"]["location[ft]"];
    m_thrust_location_x = m_thrust_location[0];
    m_thrust_location_y = m_thrust_location[1];
    m_thrust_location_z = m_thrust_location[2];
    m_T0 = data["aircraft"]["thrust"]["T0[lbf]"];
    m_T1 = data["aircraft"]["thrust"]["T1[lbf-s/ft]"];
    m_T2 = data["aircraft"]["thrust"]["T2[lbf-s^2/ft^2]"];
    m_a = data["aircraft"]["thrust"]["a"];
    get_atmospheric_properties_english(0.0, m_atm);
    m_rho0 = m_atm.density;
    //Initial
    m_airspeed = data["initial"]["airspeed[ft/s]"];
    m_altitude = data["initial"]["altitude[ft]"];
    m_heading_angle = data["initial"]["heading[deg]"];
    //Aerodynamics
    //CL
    m_CL0 = data["aerodynamics"]["CL"]["0"];
    m_CLa = data["aerodynamics"]["CL"]["alpha"];
    m_CLqbar = data["aerodynamics"]["CL"]["qbar"];
    m_CLde = data["aerodynamics"]["CL"]["de"];
    //CS
    m_CSb = data["aerodynamics"]["CS"]["beta"];
    m_CSpbar = data["aerodynamics"]["CS"]["pbar"];
    m_CSLpbar = data["aerodynamics"]["CS"]["Lpbar"];
    m_CSrbar = data["aerodynamics"]["CS"]["rbar"];
    m_CSda = data["aerodynamics"]["CS"]["da"];
    m_CSdr = data["aerodynamics"]["CS"]["dr"];
    //CD
    m_CDL0 = data["aerodynamics"]["CD"]["L0"];
    m_CDL = data["aerodynamics"]["CD"]["L"];
    m_CDL2 = data["aerodynamics"]["CD"]["L2"];
    m_CDS2 = data["aerodynamics"]["CD"]["S2"];
    m_CDqbar = data["aerodynamics"]["CD"]["qbar"];
    m_CDLqbar = data["aerodynamics"]["CD"]["Lqbar"];
    m_CDde = data["aerodynamics"]["CD"]["de"];
    m_CDLde = data["aerodynamics"]["CD"]["Lde"];
    m_CDde2 = data["aerodynamics"]["CD"]["de2"];
    //Cl
    m_Clb = data["aerodynamics"]["Cl"]["beta"];
    m_Clpbar = data["aerodynamics"]["Cl"]["pbar"];
    m_Clrbar = data["aerodynamics"]["Cl"]["rbar"];
    m_ClLrbar = data["aerodynamics"]["Cl"]["Lrbar"];
    m_Clda = data["aerodynamics"]["Cl"]["da"];
    m_Cldr = data["aerodynamics"]["Cl"]["dr"];
    //Cm
    m_Cm0 = data["aerodynamics"]["Cm"]["0"];
    m_Cma = data["aerodynamics"]["Cm"]["alpha"];
    m_Cmqbar = data["aerodynamics"]["Cm"]["qbar"];
    m_Cmde = data["aerodynamics"]["Cm"]["de"];
    //Cn
    m_Cnb = data["aerodynamics"]["Cn"]["beta"];
    m_Cnpbar = data["aerodynamics"]["Cn"]["pbar"];
    m_CnLpbar = data["aerodynamics"]["Cn"]["Lpbar"];
    m_Cnrbar = data["aerodynamics"]["Cn"]["rbar"];
    m_Cnda = data["aerodynamics"]["Cn"]["da"];
    m_CnLda = data["aerodynamics"]["Cn"]["Lda"];
    m_Cndr = data["aerodynamics"]["Cn"]["dr"];

    if (data["initial"]["type"] == "trim") {
        initalizeTrim();

    }
    else {
        initalizeState();
    }
}

void AAircraft::GetAircraftStatesUE(float& u, float& v, float& w, float& p, float& q, float& r, float& xf, float& yf, float& zf, float& e0, float& ex, float& ey, float& ez, float& Mach)
{   
    
    u = float(y[0] * 30.48);
    v = float(y[1] * 30.48);
    w = float(y[2] * 30.48);
    p = float(y[3]);
    q = float(y[4]);
    r = float(y[5]);
    xf = float(y[6] * 30.48);
    yf = float(y[7] * 30.48);
    zf = float(y[8] * 30.48);
    e0 = float(y[9]);
    ex = float(y[10]);
    ey = float(y[11]);
    ez = float(y[12]);
    get_atmospheric_properties_english(-y[8], m_atm);
    Mach = float(sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]) / m_atm.speed_of_sound);
    double V = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);

    UE_LOG(LogTemp, Warning, TEXT("states at initial = %f\t"),V);

    //UE_LOG(LogTemp, Warning, TEXT("states at initial = %f\t%f\t%f\t%f\t % f\t % f\t % f\t % f\t % f\t % f\t % f\t % f\t % f"), y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7], y[8], y[9], y[10], y[11], y[12]);
}

void AAircraft::GetAircraftControls(float& da, float& de, float& dr, float& tau)
{
    da = float(m_control[1]);
    de = float(m_control[0]);
    dr = float(m_control[2]);
    tau = float(m_control[3]);
}

void AAircraft::SetAircraftControls(float da, float de, float dr, float tau)
{
    m_control[1] = double(da);
    m_control[0] = double(de);
    m_control[2] = double(dr);
    if (tau > 1) {
        m_control[3] = 1.0;
    }
    else if (tau < 0) {
        m_control[3] = 0.0;
    }
    else {
        m_control[3] = double(tau);
    }

}

void AAircraft::TickAircraftStates(float DeltaTime)
{
    // initial state
    double ans[13];
    double rk4sol[13];
    double FM[6];

    //UE_LOG(LogTemp, Warning, TEXT("states at initial =  % f\t % f\t % f\t % f"),m_control[0], m_control[1], m_control[2], m_control[3]);

    array_copy(y, ans, 13);

    double dt = double(DeltaTime);

    

    double t0 = 0;

    rk4Solver(t0, ans, dt, 13, FM, rk4sol);

    //UE_LOG(LogTemp, Warning, TEXT("states after rk4 = %f\t%f\t%f\t%f\t % f\t % f\t % f\t % f\t % f\t % f\t % f\t % f\t % f"), ans1[0], ans1[1], ans1[2], ans1[3], ans1[4], ans1[5], ans1[6], ans1[7], ans1[8], ans1[9], ans1[10], ans1[11], ans1[12]);


    quat_norm(&rk4sol[9]);

    array_copy(rk4sol, y, 13);

}

void AAircraft::TickLatitudeLongitude(float dx, float dy, float dz, float H1)
{
}
