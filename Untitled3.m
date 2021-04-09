% B =[
% 
%     0.7073   -0.7073   -3.4956   -3.0013    3.0013    3.4956    2.1103;
%     1.1204    1.1204   -0.7919   -1.2614   -1.2614   -0.7919    0.0035;
%    -0.3309    0.3309   -0.1507   -0.3088    0.3088    0.1507   -1.2680];
% uMin =[
% 
%    -0.9599;
%    -0.9599;
%    -0.5236;
%    -0.5236;
%    -0.5236;
%    -0.5236;
%    -0.5236];
% uMax =[
% 
%     0.4363;
%     0.4363;
%     0.5236;
%     0.5236;
%     0.5236;
%     0.5236;
%     0.5236];
% yd =[0;0;2];
yd =[

   -0.0429;
   -0.0668;
   -0.0000];
   
   

uMin =[

    -0.0318;
   -0.0349;
   -0.0262;
   -0.0262;
   -0.0262;
   -0.0262;
   -0.0175];
   
   uMax =[

    0.0349;
    0.0349;
         0;
         0;
    0.0262;
    0.0262;
    0.0175];

B =[

    0.7073   -0.7073   -3.4956   -3.0013    3.0013    3.4956    2.1103;
    1.1204    1.1204   -0.7919   -1.2614   -1.2614   -0.7919    0.0035;
   -0.3309    0.3309   -0.1507   -0.3088    0.3088    0.1507   -1.2680];
[n,m] = size(B);
ep=0.1;
w=ones(m,1);
up=zeros(m,1);
A=[B           -B         -yd      zeros(n,m) zeros(n,m) zeros(n,1);
   eye(m)    zeros(m,m)  zeros(m,1) eye(m)    zeros(m,m) zeros(m,1);
   zeros(m,m) eye(m)     zeros(m,1) zeros(m,m) eye(m)    zeros(m,1);
   zeros(1,m) zeros(1,m) 1          zeros(1,m) zeros(1,m) 1];
b=[[0;0;0]-B*up;uMax-up;up-uMin;1];
c=[ep*w;ep*w;-1;zeros(m,1); zeros(m,1); 0];
% SB_LPCA(yd,B,w,up,uMin,uMax,500)
% inq = zeros((2*m+n+1),1);
%-------------------dont know--------------------------------------------
% p = revised(A,b,c,inq,'min');
% p.solve;
%---------------------ok----------------------
[x,fval,exitflag,output,lambda] = linprog(c',[],[],A,b,zeros(size(A,2),1),ones(size(A,2),1)*Inf)
%--------------------
%  u=p.x(1:m)-p.x(m+1:2*m)
%  lam=p.x(2*m+1)
 u=x(1:m)-x(m+1:2*m)
 lam=x(2*m+1)
 B*u

    
%Transform Solution Back Into control variables
% Note that x(1:m) are the + differences from the preferred control
% solution and x(m+1:m) are the - differences.

u = xout(1:m)-xout(m+1:2*m)+up;
return;
end
