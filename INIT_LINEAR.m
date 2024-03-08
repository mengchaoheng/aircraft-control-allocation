%% LOAD TRIM DATA
% ADMIRE tools used to create linear model trimmed in SSL M=0.3 at 2000 m
% load('Trim_M0p3ALT2000_LinDATA')
% ADMIRE tools used to create linear model trimmed in SSL M=0.22 at 20 m
load('Trim_M0p22ALT20_LinDATA')

%% DEFINE CONSTANTS
dt=0.01; %sim time step, seconds
d2r=pi/180; %degrees to radians conversion factor
r2d=180/pi; % radians to degrees conversion factor

%% PILOT INCEPTORS
% Fas is lateral stick force in Newtons 
% Limited to be within [-80 to 80]
Fas_v = [0 0 1 1 0 0]*40*0;
Fas_t = [0 1 1+dt 2. 2.+dt 8];
% Fes is longitudinal stick force in Newtons 
% Limited to be within [-40 to 80]
Fes_v = [0 0 1 1 -1 -1 0 0]*40;
Fes_t = [0 1 1+dt 4 4+dt 7 7+dt 8];
% Frp is rudder pedal force in Newtons 
% Limited to be within [-200 to 200]
Frp_v = [0 0 1 1 0 0]*100*0;
Frp_t = [0 1 1+dt 6. 6.+dt 8];
% TSS is steady state throttle setting in percent 
% Limited to be within [0 to 1]
TSS_v = [0 0 .1 .1 0 0]*0+u0new(10);
TSS_t = [0 1 1+dt 6. 6.+dt 8];

%% COMMANDS and REGULATORS
% Command Gradients
a1=(180*d2r)/(80^3); % Cubic roll cmd: 80 Newtons = 180 deg/s
a2=(24*d2r)/(80^3); % Cubic pitch cmd: 80 Newtons = 24 deg/s
a3=(5*d2r)/(200); % Linear yaw cmd: 200 Newtons = 5 deg
% Roll Axis Roll Rate Command
% tau_r = 0.5
Kp=2.;
% Roll Axis Bank Angle Command, not used currently
Kfp=-0.8;
Kip=0.8;
Kpp=1.6;
% Pitch Axis Pitch Rate Command
% wn=2, zeta = 0.85, Ttheta2=0
Kfq=-3.4;
Kiq=4;
Kpq=3.4;
% Yaw Axis Beta Command
% wn=6, pole/zero at 2 rad/s
Kfb=-1.71;
Kib=5.14;
Kpb=4.29;
Kr=14;

%% ON-BOARD MODEL
% For Linear DI
% Define C matrix to select control variables
% Pb, Qb, Rb
CCV=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
CCV=[CCV zeros(3,22)];
% Define CAobm to calculate the control variable derivatives
% Pdt, Qdt, Rdt
CAobm=CCV*Abare;
CAobm=CAobm(:,1:12); % select only the parts related to the 1st 12 states

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

% UsePref=1, Use the nullspace to restore towards preferred values
% UsePref=0, Do not use the nullspace to restore towards preferred values
UsePref=0; 

% CBobm only has control variables (3xm)
CBobm=CCV*Bbare;
% Define the nominal effectors as the trim settings
Unom=u0new;

global NumU Wp
NumU=16; % Number of controls
% Weighted Pseudo Inverse
W=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]); %Diagonal weighting matrix
Wp=W'*W; %Initialize Weighting matrix to unity to start
% Upref, a vector of prefered values for control solutions
Upref=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

% Ganged Controls
% Define G to make 3 virtual controls:
% d_canard = 0.5(drc+dlc)
% d_aileron = 0.5*(dloe-droe)
% d_rudder = dr
G=[ 1 0 0
    1 0 0
    0 -1 0
    0 0 0
    0 0 0
    0 1 0
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

% Define Active Control Suite, set =1 for active, =0 for inactive
% Recommend only using aerodynamic controls (1 through 7)  
UACTIVE=[
    1   % 1, Right Canard, rc, <default active>
    1   % 2, Left Canard, lc, <default active>
    1   % 3, Right Outboard Elevon, roe, <default active>
    1   % 4, Right Inboard Elevon, rie, <default active>
    1   % 5, Left Outboard Elevon, lie, <default active>
    1   % 6, Left Outboard Elevon, loe, <default active>
    1   % 7, Rudder, rud, <default active>
    0   % 8, Leading edge flaps, <default inactive>
    0   % 9, Landing Gear, <default inactive>
    0   %10, Thrust Command, <default inactive>
    0   %11, Yaw Thrust Vectoring, <default inactive>
    0   %12, Pitch Thrust Vectoring, <default inactive>
    0   %13, Forward Velocity Disturbance, <default inactive>
    0   %14, Side Velocity Disturbance, <default inactive>
    0   %15, Vertical Velocity Disturbance, <default inactive>
    0   %16, Roll Rate Disturbance, <default inactive>
];    

% Aero surface position limits are a function of Mach number.
% This data was taken from the file act_pos_lim.c provided with the ADMIRE
% simulation.  
% Note: Some of these values differ from those documented in Figure 2.2 of
% FOI-R--1624--SE.pdf 
M_vec=[0 0.5 0.8 0.9 0.95 1.4 1.5 2.5];
% Set limits based on trim condition Mach number
% Lower limits
Canard_min=interp1(M_vec,[-55 -55 -25 -15 -15 -15 -15 -15]*d2r,Mach,'next');
Elevon_min=interp1(M_vec,[-30 -30 -30 -30 -30 -30 -30 -30]*d2r,Mach,'next');
Rudder_min=interp1(M_vec,[-30 -30 -25 -20 -15 -15 -15 -10]*d2r,Mach,'next');
LEF_min=interp1(M_vec,[-10 -10 -10 -10 -10 -10 -10 -10]*d2r,Mach,'next');
% Upper limits
Canard_max=interp1(M_vec,[25 25 25 15 15 15 15 15]*d2r,Mach,'next');
Elevon_max=interp1(M_vec,[30 30 25 25 25 25 25 25]*d2r,Mach,'next');
Rudder_max=interp1(M_vec,[30 30 25 20 15 15 15 10]*d2r,Mach,'next');
LEF_max=interp1(M_vec,[30 30 30 20 20 15 10 10]*d2r,Mach,'next');
% Other Effector Position Limits
% NOTE: No position limit data available in Admire simulation for 
% disturbance effectors.  Values presented hare are place holder values.
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
% Define control limits
Umax=[Canard_max*[1 1] Elevon_max*[1 1 1 1] Rudder_max LEF_max...
    LG_max TSS_max DTY_max DTZ_max UD_max VD_max WD_max PD_max]';
Umin=[Canard_min*[1 1] Elevon_min*[1 1 1 1] Rudder_min LEF_min...
    LG_min TSS_min DTY_min DTZ_min UD_min VD_min WD_min PD_min]';


%% AIRFRAME
% For Simplified Linear Model of Aircraft
Aac=Abare(1:12,1:12); % The first 12 states
Bac=Bbare(1:12,:);
Cac=Cbare(1:12,1:12);
Dac=Dbare(1:12,:);

%% SENSORS
% SensorIDX=[1:12]; % Use "true" values of sensed variables
SensorIDX=[32:43]; % Use sensor models to estimate values of variables

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

