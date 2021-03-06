% Separate the control allocation part of INIT_NDI as a simulation of the 
% control allocation algorithm of "aircraft control allocation" to compare
% the performance of different allocation algorithms

% this is the init file
%-------------------------------------------------------------------------
% load('U1');
% U1(1:3,:)=U1(1:3,:)*0.5;
% U1(3,:)=U1(3,:)*0.5;
% load('t1');
% load('M_des');  % Pseudo control command for a certain flight,that is come from 
% M_des=zeros(100001,3);
% for i=1:100001
%     M_des(i,:)=CVdt_des(:,:,i)';
% end
%-------------------------------------------------------------------------
% aaa=randn(3,100001);
% copy from INIT_DNI,m
%-------------------------------------------------------------------------
%% LOAD TRIM DATA
% ADMIRE tools used to create linear model trimmed in SSL M=0.3 at 2000 m
% load('Trim_M0p3ALT2000_LinDATA')
% ADMIRE tools used to create linear model trimmed in SSL M=0.22 at 20 m
% load('Trim_M0p22ALT20_LinDATA')
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
% 2 - Direction Preserving==9   %无re时全局可回零，逐帧不可回零，恢复模块不但可以回零还可以使曲线圆滑
% 3 - Dirction Preserving scaled %全局局部不需要re均可回零，曲线比2光滑但re模块更佳
% 4 - Mixed Optimization, single branch
% 5 - Single Branch==8  %无re时全局可回零，逐帧不可回零，恢复模块不但可以回零还可以使曲线圆滑
% 6 - priority based on Single Branch %超幅值不能回零,解方程有问题。调节w无改善，5同。%次选
% 7 - priority based on Direction Preserving % 首选 次选
% 8 
% 9 
% 10 - priority based onDirction Preserving scaled % 首选 
% 问题，无法分解deta_m，因此优先级无法逐帧运行，因此选不考虑速度约束，恢复模块做更优解。因此5、7、10均可
LPmethod=10;

% Aero Surface Position Limits Function of Mach number
% This data was taken from the file act_pos_lim.c provided with the ADMIRE
% simulation.  
% Note: Some of these values differ from those documented in Figure 2.2 of
% FOI-R--1624--SE.pdf 


% Set limits based on trim condition Mach number
% Lower limits
% u_min=(-20)*pi/180;
u_min=umin(1);
% Upper limits
% u_max=20*pi/180;
u_max=umax(1);
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
Urlim=[250*d2r   % R canard, rad/s
    250*d2r      % L canard, rad/s
    250*d2r     % RO elevon, rad/s
    250*d2r     % RI elevon, rad/s
    250*d2r     % LI elevon, rad/s
    250*d2r     % LO elevon, rad/s
    250*d2r     % Rudder, rad/s
    210*d2r      % Leading edge flaps, rad/s
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

% % For Linear DI
% % Define C matrix to select control variables
% % Pb, Qb, Rb
% CCV=[0 0 0 1 0 0;0 0 0 0 1 0;0 0 0 0 0 1];
% CCV=[CCV zeros(3,22)];
% % CB only has control variables (3xm)
% CB=CCV*Bbare;

effector=4;
global NumU  Wp_aca
NumU=16; % Number of controls
% Weighted Pseudo Inverse
W_aca=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]); %Diagonal weighting matrix
Wp_aca=W_aca'*W_aca; %Initialize Weight matrix to unity to start
INDX=zeros(1,NumU);uMin=zeros(NumU,1);uMax=zeros(NumU,1);
CB=zeros(3,NumU);
INDX(1:effector)=ones(1,effector);
u0new=zeros(NumU,1);

uMin(1:effector)=umin;
uMax(1:effector)=umax;

% uMin(1:effector)=ones(effector,1)*(-20)*d2r;
% uMax(1:effector)=ones(effector,1)*20*d2r;

% uMax(1:effector)=ones(effector,1)*30*d2r;
% uMax(effector+1:end)=[30*d2r;LG_max; TSS_max; DTY_max; DTZ_max; UD_max; VD_max; WD_max; PD_max];
% uMin(1:effector)=ones(effector,1)*(-30)*d2r;
% uMin(effector+1:end)=[-10*d2r;LG_min; TSS_min; DTY_min; DTZ_min; UD_min; VD_min; WD_min; PD_min];

global B
CB(1:3,1:effector)=B;

% CB(1:3,1:effector)=[-0.5     0       0.5     0;
%                      0      -0.5     0       0.5;
%                      0.25    0.25    0.25    0.25];
% CB(:,1:effector) =[0.7073   -0.7073   -3.4956   -3.0013    3.0013    3.4956    2.1103;
%     1.1204    1.1204   -0.7919   -1.2614   -1.2614   -0.7919    0.0035;
%    -0.3309    0.3309   -0.1507   -0.3088    0.3088    0.1507   -1.2680];



% % Use only first 7 (aerodynamic) controls
CB2=CB(:,1:effector);
P=pinv(CB2);
% % Generate Null-Space Projection matrix
N=eye(effector,effector)-P*CB2; 
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

