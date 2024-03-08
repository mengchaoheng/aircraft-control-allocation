%% LOAD TRIM DATA
% ADMIRE tools used to create linear model trimmed in SSL M=0.3 at 2000 m
% load('Trim_M0p3ALT2000_LinDATA')
% ADMIRE tools used to create linear model trimmed in SSL M=0.22 at 20 m
load('Trim_M0p22ALT20_LinDATA')
% Clear functions with persistent variables

%% DEFINE CONSTANTS
dt=0.01; %sim time step, seconds
d2r=pi/180;
r2d=180/pi;

%% AIRCRAFT DATA
Sref=45; %Wing area, m^2
bref=10; %Wing span, m
cref=5.2; %mean chord, m
mass=9100; %kg
Ix=21000; %kg-m^2
Iy=81000; %kg-m^2
Iz=101000; %kg-m^2
Ixz=2500; %kg-m^2
GG=Ix*Iz-Ixz^2;
C1=((Iy-Iz)*Iz-Ixz^2)/GG;
C2=(Ix-Iy+Iz)*Ixz/GG;
C3=Iz/GG;
C4=Ixz/GG;
C5=(Iz-Ix)/Iy;
C6=Ixz/Iy;
C7=1/Iy;
C8=(Ix*(Ix-Iy)+Ixz^2)/GG;
C9=Ix/GG;
CC=[C3 0 C4;0 C7 0;C4 0 C9];

%% COMMANDS and REGULATORS
% Define Control Law Gains
% Command Gradients
a1=(180*d2r)/(80^3); % Cubic roll cmd: 80 Newtons = 180 deg/s
a2=(24*d2r)/(80^3); % Cubic pitch cmd: 80 Newtons = 24 deg/s
a3=(5*d2r)/(200); % Linear yaw cmd: 200 Newtons = 5 deg

%% Define mode and method constants
% Roll Mode Regulator Options
% 1 - Roll Rate Command, p
% 2 - Bank Angle Command, phi
% 3 - Crosstrack Ccommand, XTK, use with waypoint navigation
Rmode = 3;
% Pitch Mode Regulator Options
% 1 - Pitch Rate Command, q
% 2 - Flight Path Angle Command, gamma
% 3 - Altitude Command, ALT, use with waypoint navigation
Pmode = 3;

%% CONTROL ALLOCATION
% Specify Control Allocation Method
% 0 - Ganged into 3 pseudo-effectors
% 1 - Weighted Pseudo Inverse, clipped
% 2 - Weighted Pseudo Inverse, scaled
% 3 - Direct Allocation
% 4 - Cascading Generalized Inverse
% 5 - Vertex Jumping Algorithm
% 6 - Linear Programming
CAmethod=0;
% Specify which Linear Programming Algorithm
% when CAmethod is 6
% 0 - Dual Branch both 1-norm
% 1 - Dual Branch, 1-norm and inf-norm
% 2 - Direction Preserving
% 3 - Dirction Preserving scaled
% 4 - Mixed Optimization, single branch
% 5 - Single Branch
LPmethod=0;

% Aero Surface Position Limits Function of Mach number
% This data was taken from the file act_pos_lim.c provided with the ADMIRE
% simulation.  
% Note: Some of these values differ from those documented in Figure 2.2 of
% FOI-R--1624--SE.pdf 
M_vec=[0 0.5 0.8 0.9 0.95 1.4 1.5 2.5];
% Set limits based on trim condition Mach number
% Lower limits
Canard_min=[-55 -55 -25 -15 -15 -15 -15 -15]*d2r;
Elevon_min=[-30 -30 -30 -30 -30 -30 -30 -30]*d2r;
Rudder_min=[-30 -30 -25 -20 -15 -15 -15 -10]*d2r;
LEF_min=[-10 -10 -10 -10 -10 -10 -10 -10]*d2r;
% Upper limits
Canard_max=[25 25 25 15 15 15 15 15]*d2r;
Elevon_max=[30 30 25 25 25 25 25 25]*d2r;
Rudder_max=[30 30 25 20 15 15 15 10]*d2r;
LEF_max=[30 30 30 20 20 15 10 10]*d2r;
% Other Effector Position Limits
LG_min = 0; % Landing gear, 0 = Gear Up
LG_max = 1; % Landing gear, 1 = Gear Down
TSS_min = 0; % Thrust command, 0 = idle
TSS_max = 1; % Thrust command, 1 = Full Afterburner, 0.8 = max no A/B
DTY_min = -0.002; % Yaw thrust vector angle, rad (should be -25*d2r)
DTY_max = 0.002; % Yaw thrust vector angle, rad (should be 25*d2r)
DTZ_min = -0.002; % Pitch thrust vector angle, rad (should be -25*d2r)
DTZ_max = 0.002; % Pitch thrust vector angle, rad (should be -25*d2r)
UD_min = -10; % Forward velocity disturbance, m/s [No data available]
UD_max = 10; % Forward velocity disturbance, m/s [No data available]
VD_min = -10; % Side velocity disturbance, m/s [No data available]
VD_max = 10; % Side velocity disturbance, m/s [No data available]
WD_min = -10; % Vertical velocity disturbance, m/s [No data available]
WD_max = 10; % Vertical velocity disturbance, m/s [No data available]
PD_min = -10; % Roll rate disturbance, rad/s [No data available]
PD_max = 10; % Roll rate disturbance, rad/s [No data available]
% NOTE: No position limit data available in Admire simulation for 
% disturbance effectors.  Values presented hare are place holder values.

% Effector rate limits
Urlim=[50*d2r   % R canard, rad/s
    50*d2r      % L canard, rad/s
    150*d2r     % RO elevon, rad/s
    150*d2r     % RI elevon, rad/s
    150*d2r     % LI elevon, rad/s
    150*d2r     % LO elevon, rad/s
    100*d2r     % Rudder, rad/s
    10*d2r      % Leading edge flaps, rad/s
    1/dt        % Landing Gear, %/second [no rate limits in Admire]
    100/dt      % Thrust Command, %/second [no rate limits in Admire]
    1           % Yaw Thrust Vectoring, rad/s [no rate limits in Admire]
    1           % Pitch Thrust Vectoring, rad/s [no rate limits in Admire]
    1           % Forward Velocity Disturbance, m/s^2 [no rate limits in Admire]
    1           % Side Velocity Disturbance, m/s^2 [no rate limits in Admire]
    1           % Vertical Velocity Disturbance, m/s^2 [no rate limits in Admire]
    1           % Roll Rate Disturbance, rad/s^2 [no rate limits in Admire]
    ];
% NOTE: No rate limit data available in Admire simulation for effectors 8
% through 16, made up place holder values provided here. KAB

UseRL=0; % Position limited commands, UseRL=0, add rate limits UseRL=1
% Aero Actuator Dynamics
Uw=[20      % R canard
    20      % L canard
    20      % RO elevon
    20      % RI elevon
    20      % LI elevon
    20      % LO elevon
    20      % Rudder
    20];    % Leading edge flap
Mgain=diag(1./(1-exp(-Uw(1:7)*dt)));
% Delta u matrix used to perturb effectors to calculate control
% effectiveness matrix
% This is not used in the linear model, which assumes B is constant
du=[eye(7,7)*d2r eye(7,7)*(-d2r);zeros(20,14)];
% For Linear DI
% Define C matrix to select control variables
% Pb, Qb, Rb
CCV=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
CCV=[CCV zeros(3,22)];
% CB only has control variables (3xm)
CB=CCV*Bbare;
global NumU Wp
NumU=16; % Number of controls
% Weighted Pseudo Inverse
W=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]); %Diagonal weighting matrix
Wp=W'*W; %Initialize Weight matrix to unity to start

% Use only first 7 (aerodynamic) controls
CB2=CB(:,[1:7]);
P=pinv(CB2);
% Generate Null-Space Projection matrix
N=eye(7,7)-P*CB2; 
% Preferred Values
Upref=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
Upref_split=[[0 0 30 -30 -30 30 0]*d2r*.4 0 0 0 0 0 0 0 0 0]';
% Ganged Controls
% Define G to make 3 virtual controls:
% d_canard = 0.5(drc+dlc)
% d_aileron = 0.25*(dloe+dlie-drie-droe)
% d_rudder = dr
G=[ 0 1 0
    0 1 0
   -1 0 0
   -1 0 0
    1 0 0
    1 0 0
    0 0 1
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    0 0 0
    ];
% Get square B-matrix
BG=CB*G;
% Invert it
BGI=inv(BG);
% Make pseudo inverse for ganged effectors
PG=G*BGI;
% Define control limits (low speed)
Umax=[25 25 30 30 30 30 30]'*d2r;
Umin=[-55 -55 -30 -30 -30 -30 -30]'*d2r;
ulim=[Umin Umax];


% Waypoint Data for navigation
% Distances in meters
% Simulation starts at [0 0 -Alt]
% Axis is NED
Xwpt=[0 4000 4000 0 0]*8; % X distance (m), North
Ywpt=[0 0 4000 4000 0]*8; % Y distance (m), East
% Xwpt=[0 4000 4000 0 0]*1; % X distance (m), North
% Ywpt=[0 0 4000 4000 0]*1; % Y distance (m), East
% Ywpt=[0 0 -4000 -4000 0]*8; % Y distance (m), East
Wpts=[0 1 2 3 4]; % Vector of points
nwpt=5; % number of way points
WptTol=100; % Tolerance for waypoints, distance squared (m^2)
% Mission Manager Data
% Sweep the envelope
Altcmd_t=[ 0    1   100 200  300  425  450  550  600  700  800  880  980];  
Altcmd_v=[20  20   20   1500 1500 1500 1500 3000 3000 4000 4000 6000 6000];
Vcmd_t=[  0    1   100  200  300  425  450  550  600  700  800  880  980];  
Mcmd_v=[ 0.22 0.22 1.2  1.2  0.5  0.22 0.22 0.22 1.2  1.2  0.55 0.55 1.2 ];
% [t,a,p,r]=atmosisa(Altcmd_v);
a=[340 340 340 334 334 334 334 328 328 324 324 316 316]; %Speed of sound
Vcmd_v=Mcmd_v.*a;

% Bank Angle Command Vector
PHIcmd_t=[0 1 1+dt 30 30+dt 60 60+dt 90 90+dt 120 120+dt 150 150+dt 180 180+dt 210 210+dt 240 240+dt 270 270+dt 300 300+dt 330 330+dt];
PHIcmd_v=[0 0 15   15 -15   -15 0    0   30   30  -30    -30 0      0   45     45   0     0    -45   -45   0     0   -60   -60 0];

% Yaw Regulator Gains
kcb=2.5;
kib=0.2;
kpb=0.8;
kfb=-0.4;

% Attempt to decrease gain on TV
CB(:,11)=4*CB(:,11);

%% Names of signals in state-space model
%
% xdt=A*x+B*u
% y=C*x+D*u
%
% ADMIRE Tools
% Linearization has 16 inputs, 28 states and 59 outputs
%
% The variables below contain cells with the names of the various vector
% elements.
%
% The first 7 input values are the usual suite of aerodynamic control 
% effectors. 
% The 16 input values are:
Unames=[
    {'drc'}     %  1 right canard (rad)
    {'dlc'}     %  2 left canard (rad)
    {'droe'}    %  3 right outboard elevon (rad)
    {'drie'}    %  4 right inboard elevon (rad)
    {'dlie'}    %  5 left inboard elevon (rad)
    {'dloe'}    %  6 left outboard elevon (rad)
    {'dr'}      %  7 rudder (rad)
    {'dle'}     %  8 leading edge flaps (rad)
    {'ldg'}     %  9 landing gear (0-1) [0 - gear up, 1 - gear down]
    {'tss'}     % 10 thrust command (0-1)
    {'dty'}     % 11 yaw thrust vectoring (rad)
    {'dtz'}     % 12 pitch thrust vectoring (rad)
    {'u_dist'}  % 13 forward velocity disturbance (m/s)
    {'v_dist'}  % 14 side velocity disturbance (m/s)
    {'w_dist'}  % 15 vertical velocity disturbance (m/s)
    {'p_dist'}];% 16 roll rate disturbance (rad/s)
% The first 12 states are the usual ones use for aircraft dynamics. States 
% 13-28 are related to sensor dynamics. 
% The 28 states are:
Xnames=[
    {'Vt'}      %  1 velocity (m/s)
    {'alpha'}   %  2 angle of attack (rad)
    {'beta'}    %  3 sideslip angle (rad)
    {'pb'}      %  4 body axis roll rate (rad/s)
    {'qb'}      %  5 body axis pitch rate (rad/s)
    {'rb'}      %  6 body axis yaw rate (rad/s)
    {'psi'}     %  7 heading angle (rad)
    {'theta'}   %  8 pitch attitude (rad)
    {'phi'}     %  9 roll attitude (rad)
    {'x'}       % 10 X (earth parallel) position (m)
    {'y'}       % 11 Y (earth parallel) position (m)
    {'z'}       % 12 Z (earth parallel) position (m)
    {'VT1'}     % 13 VT sensor
    {'AoA1'}    % 14 alpha sensor
    {'BETA1'}   % 15 beta sensor
    {'PB1'}     % 16 pb sensor 1
    {'PB2'}     % 17 pb sensor 2
    {'QB1'}     % 18 qb sensor 1
    {'QB2'}     % 19 qb sensor 2
    {'RB1'}     % 20 rb sensor 1
    {'RB2'}     % 21 rb sensor 2
    {'THT1'}    % 22 theta sensor 1
    {'THT2'}    % 23 theta sensor 2
    {'PHI1'}    % 24 phi sensor 1
    {'PHI2'}    % 25 phi sensor 2
    {'ALT1'}    % 26 alt sensor
    {'NZ1'}     % 27 nz sensor 1
    {'NZ2'}];   % 28 nz sensor 2
% The first 31 output values are the true values from the plant. Output 
% values 32-59 are values from the sensors. In some cases, sensor dynamics
% have not been modeled and the sensed values are the same as the true 
% values.
% The 59 output values are:
Ynames=[
    {'Vt'}      %  1 velocity (m/s)
    {'alpha'}   %  2 angle of attack (rad)
    {'beta'}    %  3 sideslip angle (rad)
    {'pb'}      %  4 body axis roll rate (rad/s)
    {'qb'}      %  5 body axis pitch rate (rad/s)
    {'rb'}      %  6 body axis yaw rate (rad/s)
    {'psi'}     %  7 heading angle (rad)
    {'theta'}   %  8 pitch attitude (rad)
    {'phi'}     %  9 roll attitude (rad)
    {'x'}       % 10 X (earth parallel) position (m)
    {'y'}       % 11 Y (earth parallel) position (m)
    {'z'}       % 12 Z (earth parallel) position (m)
    {'ub'}      % 13 body axis x velocity (m/s)
    {'vb'}      % 14 body axis y velocity (m/s)
    {'wb'}      % 15 body axis z velocity (m/s)
    {'uv'}      % 16 earth parallel x axis velocity (m/s)
    {'vv'}      % 17 earth parallel y axis velocity (m/s)
    {'wv'}      % 18 earth parallel z axis velocity (m/s)
    {'nz'}      % 19 load factor of the c.g. in the z axis (g)
    {'ny'}      % 20 load factor of the c.g. in the y axis (g)
    {'Mach'}    % 21 Mach number
    {'gamma'}   % 22 flight path angle (rad)
    {'CD'}      % 23 coefficient of drag
    {'CL'}      % 24 coefficient of lift
    {'CC'}      % 25 side force coefficient
    {'Cl'}      % 26 roll moment coefficient
    {'Cm'}      % 27 pitch moment coefficient
    {'Cn'}      % 28 yaw moment coefficient
    {'Fx'}      % 29 force in the x body axis (N)
    {'Fz'}      % 30 force in the z body axis (N)
    {'My'}      % 31 moment about the y body axis (N-m)
    {'Vt_s'}    % 32 sensor velocity (m/s)
    {'alpha_s'} % 33 sensor angle of attack (rad)
    {'beta_s'}  % 34 sensor sideslip angle (rad)
    {'pb_s'}    % 35 sensor body axis roll rate (rad/s)
    {'qb_s'}    % 36 sensor body axis pitch rate (rad/s)
    {'rb_s'}    % 37 sensor body axis yaw rate (rad/s)
    {'psi_s'}   % 38 sensor heading angle (rad)
    {'theta_s'} % 39 sensor pitch attitude (rad)
    {'phi_s'}   % 40 sensor roll attitude (rad)
    {'x_s'}     % 41 sensor X (earth parallel) position (m)
    {'y_s'}     % 42 sensor Y (earth parallel) position (m)
    {'z_s'}     % 43 sensor Z (earth parallel) position (m)
    {'xdt_s'}   % 44 sensor X velocity (m/s)
    {'ydt_s'}   % 45 sensor Y velocity (m/s)
    {'zdt_s'}   % 46 sensor Z velocity (m/s)
    {'xvdt_s'}  % 47 sensor XV (m/s)
    {'yvdt_s'}  % 48 sensor YV (m/s)
    {'zvdt_s'}  % 49 sensor ZV (m/s)
    {'nz_s'}    % 50 sensor Z load factor (g)
    {'ny_s'}    % 51 sensor Y load factor (g)
    {'Mach_s'}  % 52 sensor Mach
    {'gamma_s'} % 53 sensor Flight path angle (rad)
    {'CD_s'}    % 54 sensor CD, drag coefficient
    {'CL_s'}    % 55 sensor CL, lift coefficient
    {'CY_s'}    % 56 sensor CY, side force coefficient
    {'Cl_s'}    % 57 sensor Cl, roll moment coefficient
    {'Cm_s'}    % 58 sensor Cm, pitch moment coefficient
    {'Cn_s'}];  % 59 sensor Cn, yaw moment coefficient