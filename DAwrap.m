function [u] = DAwrap(IN_MAT)
% Make single input, single output version of Direct Allocation for use in
% Simulink via the MATLAB Fcn block
% IN_MAT = [B     d
%           umin' 0
%           umax' 0
%           INDX  0]
%
% 20130921  KAB mods to try to speed up execution
% 20140524  KAB Modified to include INDX, which is used to specify active
%               effectors
global NumU
% Get sizes
[k2,m1]=size(IN_MAT);
k=k2-3;
m=m1-1;
% If matrices too small, set contols to zero and return
if k<1 || m<1 || norm(IN_MAT)<1e-16 || norm(IN_MAT(1:3,end))<1e-16
    u=zeros(NumU,1);
    return
end
% Partition input matrix into component matrices
B=IN_MAT(1:k,1:m);
v=IN_MAT(1:k,end);
umin=IN_MAT(k+1,1:m)';
umax=IN_MAT(k+2,1:m)';

INDX=IN_MAT(k+3,1:m)';

% get active effectors
B_act=B(:,INDX>0.5); 
umin_act=umin(INDX>0.5);
umax_act=umax(INDX>0.5);

[u_act,a] = da_attain(B_act,v,umin_act,umax_act);

u=zeros(NumU,1);
if(~isempty(u_act))
    u(INDX>0.5,1)=u_act;
end
end

%KAB mod of files from MC Cotting
function [u_final, mom_attain]  = da_attain (B,v,umin,umax)
% KAB 21092013
% roll=v(1);
% pitch=v(2);
% yaw=v(3);
% KAB 21092013

% CAT Control Allocation Toolbox
% from Va Tech,  Dr. Wayne Durham
% Ported by Chris Cotting, NASA DFRC
%
% Main function to perform control allocation...
%
% KAB 21092013
% cmom(1) = roll;
% cmom(2) = pitch;
% cmom(3) = yaw;
cmom(1) = v(1);
cmom(2) = v(2);
cmom(3) = v(3);
% KAB 21092013

%   cmom(1) = 0.03;
%   cmom(2) = 0.01;
%   cmom(3) = 0.02;

% add mach, alpha, etc conditions later...
%   [num_eff, B, umin, umax, num, u_1, u_2] = da_B_mat_eff_lim;
[nrow,ncol]=size(B);
num_eff=ncol;
num = 1;

for i=1:num_eff-1,
    for j=i+1:num_eff;
        u_1(num) = i;
        u_2(num) = j;
        num = num + 1;
    end
end

% Determine a non-singular 2x2 partition of the control eff. matrix
% using the facet-defining controls

allocated = 0;   % 1 = TRUE, 0 = Not yet...
limit     = 1;
Attain_mom = [];
u_allo = []; U_final = [];  mom_attain = [];
kij=[0 3 2;0 0 1];
while limit <= num -1,
    u1 = u_1(limit);
    u2 = u_2(limit);
    for i = 1:2,
        for j=i+1:3,
            A(1,1) = B(i,u1);
            A(1,2) = B(j,u1);
            A(2,1) = B(i,u2);
            A(2,2) = B(j,u2);
            A_det = det(A);
            if A_det == 0.0
                warnstring = 'Redundant controls, cannot allocate';
            elseif A_det ~= 0.0
                break
            end
        end
        if A_det ~= 0.0
            break
        end
    end
    
    % Found a non-singular partition, proceed with allocation
    
% KAB 21092013
%     if i==1 & j ==2
%         k=3;
%     elseif i == 1 & j ==3
%         k=2;
%     elseif i == 2 & j ==3
%         k = 1;
%     end
    k=kij(i,j);
% KAB 21092013
    
    %Invert the partition and get the correct row of the trans. matrix, T
    
    A_inv = inv(A);
    tr1(i) = -(A_inv(1,1)*B(k,u1) + A_inv(1,2)*B(k,u2));
    tr1(j) = -(A_inv(2,1)*B(k,u1) + A_inv(2,2)*B(k,u2));
    tr1(k) = 1.0;
    for k = 1:num_eff
        BTrow(k) = tr1(1)*B(1,k) + tr1(2)*B(2,k) + tr1(3)*B(3,k);
    end
    
    % BTrow should have 0 entries corresponding to controls u1 and u2.
    % The signs of the other entries determine whether they are a maximum
    % or minimum on the given face.  If other entries are zero, we have
    % to guess whether they are are min or max.  Introduce a new system
    % of classification where 0 == min deflection, 1 == max deflection,
    % and 2 == somewhere in between...
    
    count = 0;
    i_min(u1) = 2;
    i_min(u2) = 2;
    i_max(u1) = 2;
    i_max(u2) = 2;
    for k=1:num_eff,
        if (k ~= u1) && (k ~= u2)  % KAB 21092013

            if BTrow(k) > 0.0
                uall_min(k) = umin(k);
                uall_max(k) = umax(k);
                i_min(k) = 0;
                i_max(k) = 1;
            elseif BTrow(k) < 0.0
                uall_min(k) = umax(k);
                uall_max(k) = umin(k);
                i_min(k) = 1;
                i_max(k) = 0;
            else
                i_min(k) = 2;
                i_max(k) = 2;
                count = count +1;
                integer(count) = k;
            end
        end
    end
    
    % Now try to guess "special controls" by trial and error.  For ease
    % of calculation, assume only a "small" delta in BTrow so that the
    % control entries for u1 and u2 become positive
    
    if count > 0
        for i= 1:count
            u_min(integer(i)) = umin(integer(i));
            u_max(integer(i)) = umax(integer(i));
            i_min(integer(i)) = 0;
            i_max(integer(i)) = 1;
        end
    end
    
    % Call the sub-function MATRIX.M to determine the facet verticies.
    [mat_max] = da_matrix (num_eff, B, uall_max, u1, u2, umin, umax);
    [mat_min] = da_matrix (num_eff, B, uall_min, u1, u2, umin, umax);
    
    
    % Call the sub-function DIRALLO.M to check the current facet and
    % allocate controls if possible.  If this facet cannot be used for
    % allocation (i.e. the condition with allo_flag = 0), then check
    % opposite facet.
    
    [u_allo, flag_max, i_status, Max_mom] = da_dirallo(num_eff, B, u1, u2 ...
        ,mat_max, umin, umax, ...
        i_max, cmom);
    
    
    if flag_max == 1
        U_final = u_allo;
        Attain_mom = Max_mom;
    else
        [u_allo, flag_min, i_status, Max_mom] = da_dirallo(num_eff, B, u1, u2 ...
            ,mat_min, umin, umax, ...
            i_min, cmom);
        if flag_min == 1
            U_final = u_allo;
            Attain_mom = Max_mom;
            min_flag=1;
        end
    end
    
    limit = limit + 1;   % Increase iteration counter
    if ~isempty(Attain_mom) & ~isempty(U_final)
        mom_attain = Attain_mom;
        u_final = U_final;
        allocated = 1;
        % KAB 21092013
        break
        % KAB 21092013
    end
end
if allocated ~= 1
    u_final = [];
    mom_attain = [];
end


%   u_final
%   mom_attain
end

function [ X_mat ] = da_matrix (m ,B ,uall_X, u1, u2, umin, umax)

% CAT Control Allocation Toolbox
% from Va Tech,  Dr. Wayne Durham
% Ported by Chris Cotting, NASA DFRC
%
% Form a 3x3 matrix in moment space so that the facets determined can
% be used to allocate controls.  This matrix consists of columns
% corresponding to the facet vertex vector (referenced from the
% origin, origin), and the two facet edge vectors (referenced from the
% vertex) .


% Initialize the four verticies defining the facet and the base-2
% control vector.

Cl_mom = zeros(1,3);
Cm_mom = zeros(1,3);
Cn_mom = zeros(1,3);
%    Cx_mom = zeros(1,3);  % not needed anymore...
U      = uall_X;


% Vertex #0 with free constrols at (0,0)
U(u1) = umin(u1);
U(u2) = umin(u2);
% KAB 21092013
for i = 1:3,
    for j=1:m,
        
        Cl_mom(i) = Cl_mom(i) + B(i,j)*U(j);
    end
end
% Cl_mom=B*U';
% KAB 21092013

% Vertex #1 with free constrols at (0,1)

U(u1) = umin(u1);
U(u2) = umax(u2);
% KAB 21092013
for i = 1:3,
    for j=1:m,
        Cm_mom(i) = Cm_mom(i) + B(i,j)*U(j);
    end
end
% Cm_mom=B*U';
% KAB 21092013

% Vertex #2 with free constrols at (1,0)
U(u1) = umax(u1);
U(u2) = umin(u2);
% KAB 21092013
for i = 1:3,
    for j=1:m,
        Cn_mom(i) = Cn_mom(i) + B(i,j)*U(j);
    end
end
% Cn_mom=B*U';
% KAB 21092013

% Not needed
% Vertex #3 with free constrols at (1,1)
%  U(u1) = umax(u1);
%  U(u2) = umax(u2);
%  for i = 1:3,
%    for j=1:m,
%      Cx_mom(i) = Cx_mom(i) + B(i,j)*U(j);
%    end
%  end


% Construct matricies

for j = 1:3,
    
    
    % Construct the matrix where Column 1 is the vector from the origin to
    % the vertex 0 (both U(u1) and U(u2) are fixed at their minimums),
    % Column 2 is the vector from vertex #0 to the vertex #2 (where U(u1)
    % is the varying control), and Column 3 is the vector from vertex #0 to
    % vertex #1 (where U(u2) is the varying control)
    
    X_mat(j,1) = Cl_mom(j);
    X_mat(j,2) = Cn_mom(j) - Cl_mom(j);
    X_mat(j,3) = Cm_mom(j) - Cl_mom(j);
    
end

end

function [X_uallo, X_flag, istatus, Xmom_max] = da_dirallo(m, B, u1, u2, ...
    mat_X, umin, umax, ifac_X, ...
    X_mom)
% CAT Control Allocation Toolbox
% from Va Tech,  Dr. Wayne Durham
% Ported by Chris Cotting, NASA DFRC
%
% This subfunction is the actual direct allocation implementation.  It
% requires the control effectiveness matrix, min/max control position
% limits, the number of controls, and the previously determined facet
% geometry.  To perform allocation certain requirements must be met.
% These are presented below.

X_uallo = [];
Xmom_max = [];

X_flag = 0;        % Allocation flag 0,1   1 == TRUE
istatus = 0;
mat_det = det(mat_X);
mat_inv = inv(mat_X);

% This section allocates the controls provided that the determinant of
% the matrix containging the facet geometry (mat_X) is not zero.

% if mat_det ~= 0.0    % desired condition
if abs(mat_det) > eps    % desired condition KAB Mod Never hitting else condition
    c1 = 0.0;
    c2 = 0.0;
    c3 = 0.0;
    for i = 1:3;
        c1 = c1 + mat_inv(1,i)*X_mom(i);
        c2 = c2 + mat_inv(2,i)*X_mom(i);
        c3 = c3 + mat_inv(3,i)*X_mom(i);
    end
    
    %Test the value of c1.  If it is zero, then allocation cannot
    %continue, if it is greater than zero, allocation proceeds, if it is
    %less than or equal to zero then there is negative saturation.
    
    if c1 == 0.0
        istatus = -1;   %Moment is parallet to the facet
        return
    end
    
    if c1 > 0.0
        c2 = c2/c1;
        c3 = c3/c1;
        if c1 > 1.0
            c1 = 1.0;
        end
%         if (c2 >= 0.0) & (c3 >= 0.0) & (c2 <=1.0) & (c3 <= 1.0)
% KAB Mod: add eps, because algorithm not finding solution if moment is
% pointed add edge (==1 or ==0 were off by a small precision error)
        if (c2 >= 0.0-eps) && (c3 >= 0.0-eps) && (c2 <=1.0+eps) && (c3 <= 1.0+eps) % KAB 21092013
            istatus = 0;      % Everything is good
            [X_uallo, Xmom_max] = da_allocate (m, B, u1, u2, ifac_X, c1, ...
                c2, c3, umin, umax);
            X_flag = 1;
            allo1=1  ;% KILL ME!!
        else
            istatus = 6;      % Everything is ok, just have the wrong facet
        end
    else                  % if c1 <= 0.0
        istatus = 7;        % negative saturation (a bad thing...)
    end
    
    %  This section tests to see if the origin is on the boundary, or if
    %  the edges in question are parallel.  If the origin is on the
    %  boundary, the moment may be as well and allocation is still
    %  possible.  If the edges are parallel, allocation of hte defining
    %  controls is possible.
    
% elseif mat_det == 0
else %KAB Mod
    A2(1,1) = mat_X(1,2);
    A2(1,2) = mat_X(1,3);
    A2(2,1) = mat_X(2,2);
    A2(2,2) = mat_X(2,3);
    A2_det = det(A2);
    if A2_det == 0          % Co-linear edges, cannot allocate
        return
    end
    
    % Found a non-singular 2x2 partition of a geometry matrix, proceed
    % with tests
    
%     i = 0;
%     j = 1;
%     k = 2;
% KAB I think the above was from C-code, Matlab does not index 0-2
    i = 1;
    j = 2;
    k = 3;
    
    A2_inv = inv(A2);
    
    % Determine if the defining vertex (ie column 1 of mat_X) is in the
    % place of the facet being studied.
    
    c2 = -(A2_inv(1,1)*mat_X(i,1) + A2_inv(1,2)*mat_X(j,1));
    c3 = -(A2_inv(2,1)*mat_X(i,1) + A2_inv(2,2)*mat_X(j,1));
    c1 = c2*mat_X(k,2) + c3*mat_X(k,3);
    
    if mat_X(k,1) ~= c1
        istatus = 1;            % Singular Origin not on boundary
        return
    end
    
    c2 = -(A2_inv(1,1)*X_mom(i) + A2_inv(1,2)*X_mom(j));
    c3 = -(A2_inv(2,1)*X_mom(i) + A2_inv(2,2)*X_mom(j));
    c1 = c2*mat_X(k,2) + c3*mat_X(k,3);
    
    if X_mom(k) ~= c1
        status = 3;             % Origin on boundary, moment is not
        return
    end
    
    % Now determine if the origin is on the facet or just in the plane of
    % the facet.
    
    if (c2 >= 0.0) & (c3 >= 0.0) & (c2 <= 1.0) & (c3 <= 1.0)
        % Origin and moment on boundary, can allocate
        [X_uallo, Xmom_max] = da_allocate (m, B, u1, u2, ifac_X, c1, ...
            c2, c3, umin, umax);
        X_flag = 1;
        allo2=1  % KILL ME!!
    else
        return                  % Origin on boundary, moment is not
    end
    
end

end
 
 function [U_vec, M_vec] = da_allocate (m, B, u1, u2, i_factor, var1, ...
     var2, var3, umin, umax)
 
 % CAT Control Allocation Toolbox
 % from Va Tech,  Dr. Wayne Durham
 % Ported by Chris Cotting, NASA DFRC
 %
 % This subfunction proceeds only if the control allocation is
 % determined to be possible.  The purpose of this routine is to find
 % the allocation control vector X_uallo and the corresponding maximum
 % moment vector, Xmom_max.
 
% KAB 21092013
%  for i = 1:3,
%      for j = 1:m,            % Initialize maximum moment vector
%          M_vec(i) = 0.0;
%      end
%  end
M_vec=[0 0 0];
% KAB 21092013
 
 for i = 1:m,
     if (i ~= u1) && ( i ~= u2)   % Effectively skips u1 and u2 % KAB 21092013
         if i_factor(i) == 0
             U_vec(i) = umin(i) *var1;
         elseif i_factor(i) == 1
             U_vec(i) = umax(i) * var1;
         end
     end
 end
 U_vec(u1) = var1*(umin(u1) + var2*(umax(u1) - umin(u1)));
 U_vec(u2) = var1*(umin(u2) + var3*(umax(u2) - umin(u2)));
 
 for i = 1:3,
     for j = 1:m,
         M_vec(i) = M_vec(i) + B(i,j)*U_vec(j);
     end
 end
 
 end
  
      
      