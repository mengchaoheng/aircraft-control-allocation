clear all;
close all;
global NumU Wp
NumU=16; % Number of controls
% Weighted Pseudo Inverse
W=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]); %Diagonal weighting matrix
Wp=W'*W; %Initialize Weight matrix to unity to start
B=zeros(3,NumU);uMin=zeros(NumU,1);uMax=zeros(NumU,1);INDX=zeros(1,NumU);
effector=4;
n=3;
B(:,1:effector)=[-0.50004   0       0.5001   0;
    0  -0.5005    0       0.50003;
    0.25005   0.2502   0.251   0.2503];
uMin(1:effector)=ones(effector,1)*(-20)*pi/180;
uMax(1:effector)=ones(effector,1)*20*pi/180;

% B(:,1:effector) =[0.7073   -0.7073   -3.4956   -3.0013    3.0013    3.4956    2.1103;
%     1.1204    1.1204   -0.7919   -1.2614   -1.2614   -0.7919    0.0035;
%    -0.3309    0.3309   -0.1507   -0.3088    0.3088    0.1507   -1.2680];
% uMin(1:effector) =[-0.9599;-0.9599;-0.5236;-0.5236;-0.5236;-0.5236;-0.5236];
% uMax(1:effector) =[0.4363;0.4363;0.5236;0.5236;0.5236;0.5236;0.5236];

INDX(1:effector)=ones(1,effector);
% IN_MAT = [B     d ye;
%           umin' 0 0;
%           umax' 0 0;
%           INDX  LPmethod 0];



LPmethod=8;% just 2 3 7 ok
yd=[0;0;0]; 
ye=[0;0;0];
IN_MAT = [B yd ye;uMin' 0 0;uMax' 0 0;INDX  LPmethod 0];
% [u] = LPwrap(IN_MAT);      
% u(INDX)      
% B(:,INDX>0.5)*u(INDX)        
%===================================·ùÖµ²âÊÔ==================================================
N=10;
x=zeros(effector,(N+1)^2);
u=zeros(effector,1);
[X,Y,Z] = sphere(N);
load('M_des');
% t=0:0.01:100;
% X=0.04*sin(0.01*pi*t);
% Y=0.04*cos(0.02*pi*t);
% Z=0*sin(pi*t);
x1=zeros(effector,(N+1)^2);
u1=zeros(effector,1);
% x2=zeros(effector,(N+1)^2);
% u2=zeros(effector,1);

for i=1:(N+1)^2%length(M_des(1:1000,1))%%length(X)
v=0.0000001*[X(i);Y(i);Z(i)];% ÐéÄâÖ¸ÁîM_des(i,:)'%
IN_MAT(1:3,end-1)=v;

%=====================================
% u=pinv(B)*v;
% x(:,i)=Constrain(u,uMin,uMax);

[u1] = LPwrap(IN_MAT);
x1(:,i)=u1(INDX>0.5);%Constrain(u1(INDX>0.5),uMin,uMax);

% u1=DPscaled_LPCA(v,B,uMin,uMax,100);
% x1(:,i)=Constrain(u1,uMin,uMax);

% u2=SBprio_LPCA(v,ye,B,ep*ones(4,1),zeros(4,1),uMin,uMax,5e2);
% x2(:,i)=Constrain(u2,uMin,uMax);
end
U=B(:,INDX>0.5)*x;
U1=B(:,INDX>0.5)*x1;
% U2=B*x2;
% V=(yd+ye);
figure(1),
% plot3(U(1,:),U(2,:),U(3,:),'b*');
% hold on;
plot3(U1(1,:),U1(2,:),U1(3,:),'r*');
% hold on;
% plot3(ye(1,1),ye(2,1),ye(3,1),'b*');
% hold on;
% plot3(V(1,1),V(2,1),V(3,1),'g*');
% hold on;
% plot3(u2(1,:),u2(2,:),u2(3,:),'effector*');
% 
% hold on;
% plot3(u4(1,:),u4(2,:),u4(3,:),'g>');
