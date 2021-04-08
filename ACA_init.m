% Separate the control allocation part of INIT_NDI as a simulation of the 
% control allocation algorithm of "aircraft control allocation" to compare
% the performance of different allocation algorithms

% this is the init file
clear all;
%-------------------------------------------------------------------------
load('U1');
% U1(1:3,:)=U1(1:3,:)*0.5;
% U1(3,:)=U1(3,:)*0.5;
load('t1');
load('M_des');  % Pseudo control command for a certain flight,that is come from 
% M_des=zeros(100001,3);
% for i=1:100001
%     M_des(i,:)=CVdt_des(:,:,i)';
% end
%-------------------------------------------------------------------------
aaa=randn(3,100001);
% copy from INIT_DNI,m
%-------------------------------------------------------------------------
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
%% CONTROL ALLOCATION
% Specify Control Allocation Method
% 0 - Ganged into 3 pseudo-effectors
% 1 - Weighted Pseudo Inverse, clipped
% 2 - Weighted Pseudo Inverse, scaled
% 3 - Direct Allocation
% 4 - Cascading Generalized Inverse
% 5 - Vertex Jumping Algorithm
% 6 - Linear Programming
CAmethod=6;
% Specify which Linear Programming Algorithm
% when CAmethod is 6
% 0 - Dual Branch both 1-norm
% 1 - Dual Branch, 1-norm and inf-norm
% 2 - Direction Preserving
% 3 - Dirction Preserving scaled, maybe non-restore
% 4 - Mixed Optimization, single branch
% 5 - Single Branch
% 6 - priority based on Single Branch
% 7 - priority based on Direction Preserving
LPmethod=6;

% Aero Surface Position Limits Function of Mach number
% This data was taken from the file act_pos_lim.c provided with the ADMIRE
% simulation.  
% Note: Some of these values differ from those documented in Figure 2.2 of
% FOI-R--1624--SE.pdf 


% Set limits based on trim condition Mach number
% Lower limits
u_min=(-20)*pi/180;

% Upper limits
u_max=20*pi/180;

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
Urlim=[200*d2r   % R canard, rad/s
    200*d2r      % L canard, rad/s
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
% Urlim=[350*d2r   % R canard, rad/s
%     350*d2r      % L canard, rad/s
%     350*d2r     % RO elevon, rad/s
%     350*d2r     % RI elevon, rad/s
%     1     % LI elevon, rad/s
%     1     % LO elevon, rad/s
%     1     % Rudder, rad/s
%     1      % Leading edge flaps, rad/s
%     1        % Landing Gear, %/second [no rate limits in Admire]
%     1      % Thrust Command, %/second [no rate limits in Admire]
%     1           % Yaw Thrust Vectoring, rad/s [no rate limits in Admire]
%     1           % Pitch Thrust Vectoring, rad/s [no rate limits in Admire]
%     1           % Forward Velocity Disturbance, m/s^2 [no rate limits in Admire]
%     1           % Side Velocity Disturbance, m/s^2 [no rate limits in Admire]
%     1           % Vertical Velocity Disturbance, m/s^2 [no rate limits in Admire]
%     1           % Roll Rate Disturbance, rad/s^2 [no rate limits in Admire]
%     ];
% NOTE: No rate limit data available in Admire simulation for effectors 8
% through 16, made up place holder values provided here. KAB

UseRL=0; % Position limited commands, UseRL=0, add rate limits UseRL=1


% For Linear DI
% Define C matrix to select control variables
% Pb, Qb, Rb
CCV=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
CCV=[CCV zeros(3,22)];
% CB only has control variables (3xm)
CB=CCV*Bbare;
% CB=zeros(3,16);
% CB(1:3,1:4)=[-0.5   0       0.5   0;
%     0  -0.5    0       0.5;
%     0.25   0.25   0.25   0.25];
% CB(1:3,1:7) =[1.9808   -1.9805  -20.0502  -18.0003   17.9996   20.0471    7.7784;
%     4.7436    4.7435   -5.0333   -8.3523   -8.3526   -5.0336         0;
%    -1.3682    1.3681   -0.9695   -2.4887    2.4886    0.9693   -7.2922];
global NumU Wp
NumU=16; % Number of controls
% Weighted Pseudo Inverse
W=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]); %Diagonal weighting matrix
Wp=W'*W; %Initialize Weight matrix to unity to start

% % Use only first 7 (aerodynamic) controls
CB2=CB(:,[1:7]);
P=pinv(CB2);
% % Generate Null-Space Projection matrix
N=eye(7,7)-P*CB2; 
% Preferred Values
Upref=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
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
BGI=pinv(BG);
% Make pseudo inverse for ganged effectors
PG=G*BGI;

