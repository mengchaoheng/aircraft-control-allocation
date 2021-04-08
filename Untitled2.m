% ===================================∑˘÷µ≤‚ ‘==================================================
yd =[

   -4.3016;
   -0.2706;
   -0.2706];
B =[

    0.7073   -0.7073   -3.4956   -3.0013    3.0013    3.4956    2.1103;
    1.1204    1.1204   -0.7919   -1.2614   -1.2614   -0.7919    0.0035;
   -0.3309    0.3309   -0.1507   -0.3088    0.3088    0.1507   -1.2680];
uMin =[

   -0.9599;
   -0.9599;
   -0.5236;
   -0.5236;
   -0.5236;
   -0.5236;
   -0.5236];
uMax =[

    0.4363;
    0.4363;
    0.5236;
    0.5236;
    0.5236;
    0.5236;
    0.5236];
[m,k] = size(B);
ep=0.1;
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
x2=zeros(k,(N+1)^2);
u2=zeros(k,1);

for i=1:(N+1)^2%%length(X)
v=10*[X(i);Y(i);Z(i)];% –Èƒ‚÷∏¡Ó
%=====================================
% u=pinv(B)*v;
% x(:,i)=Constrain(u,umin,umax);

% u1=DP_LPCA(v,B,umin,umax,100);
% x1(:,i)=Constrain(u,umin,umax);

u=DPscaled_LPCA(v,B,uMin,uMax,500);
x(:,i)=Constrain(u,uMin,uMax);
% u2=SBprio_LPCA(v,ye,B,ep*ones(4,1),zeros(4,1),umin,umax,5e2);
% x2(:,i)=Constrain(u2,umin,umax);
end
U=B*x;
U1=B*x1;
U2=B*x2;

figure(1),
plot3(U(1,:),U(2,:),U(3,:),'b*');
hold on;
% plot3(U1(1,:),U1(2,:),U1(3,:),'r*');
% hold on;
% plot3(ye(1,1),ye(2,1),ye(3,1),'b*');
% hold on;
plot3(yd(1,1),yd(2,1),yd(3,1),'g*');
% hold on;
% plot3(U2(1,:),U2(2,:),U2(3,:),'k*');