function [ UK ] = getUperp( IN_MAT )
%[ UK ] = getUperp( IN_MAT )
% IN_MAT = [B
%           deltaU+Unom+..  %使用全局的u，为了在B的零空间找到范数最小u
%           INDX]
%   calculates the the component of the direction towards the preferred
%   solution which lies in the null-space of B
%   B must be wider than it is tall (more columns than rows)
% 20150823  KAB Created file to implement Durham's min-norm restoring
% 20160301  KAB Added option to use weighting matrix

% get u_perpendicular, in 7.4.4.1
global NumU Wp 
[n1,m1]=size(IN_MAT);
INDX=IN_MAT(end,:)>0.5;
Bx=IN_MAT(1:n1-1,INDX);
[n,m]=size(Bx);
UK=zeros(NumU,1);
% If weighting matrix is identity
if norm(Wp(INDX>0.5,INDX>0.5)-eye(sum(INDX>0.5)))<eps
    % Minimum Norm restoring with equal weighting on effectors, min u'*u
    VV=zeros(n,1);
    VV(n,1)=-2;
    Uperp=pinv(Bx)*VV;
    if all(abs(Uperp)<=eps)
        Kopt=0; %防止除0
    else
        Kopt=2/(Uperp'*Uperp);
    end
    if Kopt>1
        UK(INDX,1)=0*Uperp; %消除抖动
    else
        UK(INDX,1)=Kopt*Uperp;
    end
else
    % Minimum Norm restoring with weighting matrix, min u'*Wp*u
    B=Bx(1:n-1,:);
    Wpp=Wp(INDX>0.5,INDX>0.5);
    P=Wpp'\B'/(B*(Wpp'\B'));
    N=eye(m,m)-P*B;
    un=-N*Bx(n,:)';
    UK(INDX,1)=un;
end

end

