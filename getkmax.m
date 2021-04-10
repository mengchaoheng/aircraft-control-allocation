function [ K ] = getkmax( IN_MAT )
%[ UK ] = getUperp( IN_MAT )
% IN_MAT = [MaxValues
%           MinValues
%           INDX
%           u_kopt]
%   calculates the the component of the direction towards the preferred
%   solution which lies in the null-space of B
%   B must be wider than it is tall (more columns than rows)
% 20150823  KAB Created file to implement Durham's min-norm restoring
% 20160301  KAB Added option to use weighting matrix

% LimitScale, get K_max, in 7.4.4.1

INDX=IN_MAT(3,:)>0.5;
u_kopt=IN_MAT(4,INDX)
MaxValues=IN_MAT(1,INDX)
MinValues=IN_MAT(2,INDX)
% if any(abs(MaxValues)<eps) || any(abs(MinValues)<eps)
%     K=0;
%     return;
% end
%---------------------test 0-----------
if any(MaxValues<0) 
    MaxValues(MaxValues<0)=eps;
end
if any(MinValues>0) 
    MinValues(MinValues>0)=-eps;
end
% u_abs=abs(u_kopt);
u_abs=[u_kopt(u_kopt>0) -u_kopt(u_kopt<0)]
Max=[MaxValues(u_kopt>0) -MinValues(u_kopt<0)]


K_temp=Max(u_abs>eps)./u_abs(u_abs>eps);
K=min([1 K_temp]);
% for i=1:size(INDX)
%     if  u_abs(i)>eps
%         K_temp=Max(i)/u_abs(i);
%         if K_temp<K
%             K=K_temp;
%         end
%     end
% end
end

