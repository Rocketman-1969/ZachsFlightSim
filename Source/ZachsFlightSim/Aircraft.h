// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "foster_helper.h"
#include "CoreMinimal.h"
#include "GameFramework/Pawn.h"
#include "Aircraft.generated.h"

UCLASS()
class ZACHSFLIGHTSIM_API AAircraft : public APawn
{
	GENERATED_BODY()

public:
	// Sets default values for this pawn's properties
	AAircraft();

protected:
    // Called when the game starts or when spawned
    virtual void BeginPlay() override;

    // declare variables
    Atmosphere m_atm;
    json data;

    //simulation
    double m_time_step;
    double m_total_time;
    double m_rho0;

    //aircraft
    double m_wing_area;
    double m_wing_span;
    double m_wing_chord;
    double m_weight;
    double m_Ixx;
    double m_Iyy;
    double m_Izz;
    double m_Ixy;
    double m_Ixz;
    double m_Iyz;
    double m_hx;
    double m_hy;
    double m_hz;
    double m_CG_shift_x;
    double m_CG_shift_y;
    double m_CG_shift_z;

    //thrust
    double m_thrust_location_x;
    double m_thrust_location_y;
    double m_thrust_location_z;
    double m_thrust_direction_x;
    double m_thrust_direction_y;
    double m_thrust_direction_z;
    double m_T0;
    double m_T1;
    double m_T2;
    double m_a;

    //initial
    double m_airspeed;
    double m_altitude;
    double m_elevation_angle;
    double m_bank_angle;
    double m_alpha;
    double m_beta;
    double m_p;
    double m_q;
    double m_r;
    double m_heading_angle;
    double m_ailoron;
    double m_elevator;
    double m_rudder;
    double m_throttle;

    // aerodynamics
    double m_CL0;
    double m_CLa;
    double m_CLqbar;
    double m_CLde;
    double m_CSb;
    double m_CSpbar;
    double m_CSLpbar;
    double m_CSrbar;
    double m_CSda;
    double m_CSdr;
    double m_CDL0;
    double m_CDL;
    double m_CDL2;
    double m_CDS2;
    double m_CDqbar;
    double m_CDLqbar;
    double m_CDde;
    double m_CDLde;
    double m_CDde2;
    double m_Clb;
    double m_Clpbar;
    double m_Clrbar;
    double m_ClLrbar;
    double m_Clda;
    double m_Cldr;
    double m_Cm0;
    double m_Cma;
    double m_Cmqbar;
    double m_Cmde;
    double m_Cnb;
    double m_Cnpbar;
    double m_CnLpbar;
    double m_Cnrbar;
    double m_Cnda;
    double m_CnLda;
    double m_Cndr;

    //Forces and Moments
    double m_control[6];
    double m_Imatrix[3][3];
    double m_hmatrix[3][3];
    double m_ImatrixInv[3][3];

    //RK4 variables
    int const m_size = 13;
    double m_k1[13];
    double m_k2[13];
    double m_k3[13];
    double m_k4[13];

    //Forces moments variables 
    //double m_FM[6];

    //initalize sim
    double y[13];


    void rk4Solver(double t, double* y0, double dt, int size, double* FM, double* ans);
    void PhysicsModel(double t, double* y0, double* FM, double* ans);
    void aerodynamicModel(double* y0, double* FM);
    void initalizeState();
    void initalizeTrim();
    void GetR(double* G, double* y0, double g, double elevation_angle, double bank_angle, double beta, bool side_slip, double* FM, double* R);

    UFUNCTION(BlueprintCallable, Category = "FlightSim")
    void InitalizeAircraftFromJSON(FString ConfigFileName);

    UFUNCTION(BlueprintCallable, Category = "FlightSim")
    void GetAircraftStatesUE(float& u, float& v, float& w, float& p, float& q, float& r, float& xf, float& yf, float& zf, float& e0, float& ex, float& ey, float& ez, float& Mach);

    UFUNCTION(BlueprintCallable, Category = "FlightSim")
    void GetAircraftControls(float& da, float& de, float& dr, float& tau);

    UFUNCTION(BlueprintCallable, Category = "FlightSim")
    void SetAircraftControls(float da, float de, float dr, float tau);

    UFUNCTION(BlueprintCallable, Category = "FlightSim")
    void TickAircraftStates(float DeltaTime);

    void TickLatitudeLongitude(float dx, float dy, float dz, float H1);

public:	
	// Called every frame
	virtual void Tick(float DeltaTime) override;

	// Called to bind functionality to input
	virtual void SetupPlayerInputComponent(class UInputComponent* PlayerInputComponent) override;

};
