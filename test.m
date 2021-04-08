clear all;
close all;
global NumU Wp
NumU=16; % Number of controls
% Weighted Pseudo Inverse
W=diag([1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]); %Diagonal weighting matrix
Wp=W'*W; %Initialize Weight matrix to unity to start
B=zeros(3,NumU);uMin=zeros(NumU,1);uMax=zeros(NumU,1);INDX=false(1,16);
k=7;
n=3;
% B=[-0.5   0       0.5   0;
%     0  -0.5    0       0.5;
%     0.25   0.25   0.25   0.25];
% uMin=[1;1;1;1]*(-20)*pi/180;
% uMax=[1;1;1;1]*20*pi/180;

B(:,1:k) =[0.7073   -0.7073   -3.4956   -3.0013    3.0013    3.4956    2.1103;
    1.1204    1.1204   -0.7919   -1.2614   -1.2614   -0.7919    0.0035;
   -0.3309    0.3309   -0.1507   -0.3088    0.3088    0.1507   -1.2680];
uMin(1:k) =[-0.9599;-0.9599;-0.5236;-0.5236;-0.5236;-0.5236;-0.5236];
uMax(1:k) =[0.4363;0.4363;0.5236;0.5236;0.5236;0.5236;0.5236];

INDX(1:k)=[true true true true true true true];

% IN_MAT = [B     d;
%           umin' 0;
%           umax' 0;
%           INDX  LPmethod];



LPmethod=7;% just 2 3 7 ok
yd=[1;0;0]; 
IN_MAT = [B yd;uMin' 0;uMax' 0;INDX  LPmethod];
% [u] = LPwrap(IN_MAT);      
% u(INDX)      
% B(:,INDX>0.5)*u(INDX)        
%===================================·ùÖµ²âÊÔ==================================================
N=100;
x=zeros(k,(N+1)^2);
u=zeros(k,1);
[X,Y,Z] = sphere(N);
load('M_des');
% t=0:0.01:1;
% X=0.3*sin(pi*t);
% Y=0.3*cos(pi*t);
% Z=0*sin(pi*t);
x1=zeros(k,(N+1)^2);
u1=zeros(k,1);
% x2=zeros(k,(N+1)^2);
% u2=zeros(k,1);

for i=1:(N+1)^2%%length(X)%length(M_des(1:1000,1))%
v=15*[X(i);Y(i);Z(i)];% ÐéÄâÖ¸ÁîM_des(i,:)'%
IN_MAT(1:n,end)=v;

%=====================================
% u=pinv(B)*v;
% x(:,i)=Constrain(u,uMin,uMax);

[u1] = LPwrap(IN_MAT);
x1(:,i)=Constrain(u1(INDX),uMin,uMax);

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
% plot3(u2(1,:),u2(2,:),u2(3,:),'k*');
% 
% hold on;
% plot3(u4(1,:),u4(2,:),u4(3,:),'g>');
