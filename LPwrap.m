function [u] = LPwrap(IN_MAT)
% Make single input, single output version of Linear Programming for use in
% Simulink via the MATLAB Fcn block
% IN_MAT = [B     d ye
%           umin' 0 0
%           umax' 0 0
%           INDX  LPmethod 0]
%
% 20140905  Created version to use Roger Beck's DB_LPCA program
% 20151206  Updated with Roger Beck's latest version of code and added 
%           LPmethod to select the various Linear Programming algorithms 

global NumU
% Get sizes
[k2,m1]=size(IN_MAT);
k=k2-3;
m=m1-2;
% If matrices too small, set contols to zero and return
if k<1 || m<1 || norm(IN_MAT)<1e-16
    u=zeros(NumU,1);
    return
end
% Partition input matrix into component matrices
B=IN_MAT(1:k,1:m);
v=IN_MAT(1:k,end-1);
umin=IN_MAT(k+1,1:m)';
umax=IN_MAT(k+2,1:m)';
LPmethod=IN_MAT(end,end-1);
ye=IN_MAT(1:k,end);
% LPmethod should be an integer between 0 and 5
% 0 = DB_LPCA
        %   Dual Branch Control Allocation - Linear Program
        %      Objective Error Minimization Branch (1-norm)
        %      Control Error Minimization (1-norm)
% 1 = DBinf_LPCA
        %   Dual Branch Control Allocation - Linear Program
        %      Objective Error Minimization (1-norm)
        %      Control Error Minimization (inf-norm)
% 2 = DP_LPCA
        % Direction Preserving Control Allocation Linear Program    
% 3 =DPscaled_LPCA
        % Direction Preserving Control Allocation Linear Program
        %     Reduced formulation (Solution Scaled from Boundary)
% 4 = MO_LPCA
        % Mixed Optimization (Single Branch) Control Allocation Linear Program
        %    Objective Error Minimizing
        %    Control Error minimizing
% 5 = SB_LPCA
        % Single Branch Control Allocation Linear Program
        %    Direction Preserving
        %    Control Error minimizing
% If LPmethod is not one of the allowed options, set it to zero
if sum(LPmethod==[0 1 2 3 4 5 6 7])~=1
    LPmethod=0;
end

INDX=IN_MAT(k+3,1:m)';
m_act=sum(INDX);
% get active effectors
B_act=B(:,INDX>0.5); 
umin_act=umin(INDX>0.5);
umax_act=umax(INDX>0.5);

%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          wd [n,1]  = Weighting on Objective error
%          up [m,1]  = Preferred control solution
%          wu [m,1]  = Weighting on control error
%          emax[n,1] = Upper Bound for objective error
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%          lam       = Control Error Weighting (single parameter)
%          eMax[n,1] = Maximum Objective error
%          w [m,1]   = Control Error Weighting

% change variable names to be consistent with Roger's documentation
yd=v;
B=B_act;
uMin=umin_act;
uMax=umax_act;
% Some variables not specifically defined for generic case
% Values below are suggestions for now
% Users may add code here to specify alternative values
wd=ones(k,1);
up=zeros(m_act,1);
wu=0.1*ones(m_act,1);
emax=ones(k,1)*2e1;
itlim=50;
lam=0.1;
eMax=emax;
w=wu;
switch LPmethod
    case 0
        [u_act, feas, errout,itlim] = DB_LPCA(yd+ye,B,wd,up,wu,emax,...
            uMin,uMax,itlim);
        %   Dual Branch Control Allocation - Linear Program
        %      Objective Error Minimization Branch (1-norm)
        %      Control Error Minimization (1-norm)
    case 1
        [u_act, feas, errout,itlim] = DBinf_LPCA(yd+ye,B,wd,up,wu,emax,...
            uMin,uMax,itlim);
        %   Dual Branch Control Allocation - Linear Program
        %      Objective Error Minimization (1-norm)
        %      Control Error Minimization (inf-norm)
    case 2
        [u_act, errout] = DP_LPCA(yd+ye,B,uMin,uMax,itlim);
        % Direction Preserving Control Allocation Linear Program    
    case 3
        [u_act,itlim,errout] = DPscaled_LPCA(yd+ye,B,uMin,uMax,itlim);
        % Direction Preserving Control Allocation Linear Program
        %     Reduced formulation (Solution Scaled from Boundary)
    case 4
        [u_act,errout] = MO_LPCA(yd+ye,B,up,lam, eMax,uMin,uMax,itlim);
        % Mixed Optimization (Single Branch) Control Allocation Linear Program
        %    Objective Error Minimizing
        %    Control Error minimizing
    case 5
        [u_act,errout] = SB_LPCA(yd+ye,B,w,up,uMin,uMax,itlim);
        % Single Branch Control Allocation Linear Program
        %    Direction Preservingyd
        %    Control Error minimizing
    case 6
        [u_act,errout] = SBprio_LPCA(yd,ye,B,w,up,uMin,uMax,itlim);
        % Single Branch Control Allocation Linear Program
        %    Direction Preserving
        %    Control Error minimizing
    case 7
        [u_act,errout] = DPprio_LPCA(yd,ye,B,uMin,uMax,itlim);
        % Single Branch Control Allocation Linear Program
        %    Direction Preserving
        %    Control Error minimizing
end
u=zeros(NumU,1);
u(INDX>0.5,1)=u_act;

end

function [u, feas, errout,itlim] = DB_LPCA(yd,B,wd,up,wu,emax,uMin,uMax,itlim) % note
%   Dual Branch Control Allocation - Linear Program
%      Objective Error Minimization Branch (1-norm)
%      Control Error Minimization (1-norm)
%
% [u, feas, errout,itlim] = DB_LPCA(yd,B,wd,up,wu,emax,uMin,uMax,itlim);
%
%    Uses a Bounded Revised Simplex solver to solve two linear Programming problems
%  The first branch, feasibility, seeks to minimize the weighted 1-norm of the
%    objective error
%    min J= | diag(wd)*(B*u - yd) |_1  s.t. umin <= u <=umax
%
%   If J = 0 then yd is on the interior of the AMS and a second program seeks to minimize
%   weighted error with a "preferred" control solution.
%    min |diag(wu)*(u-up)|_1  s.t. umin <= u <= umax and Bu=yd
%
%   (In the text, the overall structure is presented in A.5.1)
%
%   (See Buffington, J. et al. "On-Line System Identification for Aircraft with Distributed
%         Control Effectors", Int. J. Robust Nonlinear Control 9, 1033-1049 (1999)).
%        
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          wd [n,1]  = Weighting on Objective error
%          up [m,1]  = Preferred control solution
%          wu [m,1]  = Weighting on control error
%          emax[n,1] = Upper Bound for objective error
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         feas       = Feasible flag
%                        0 = yd is unachievable
%                        1 = Bu=yd
%                        2 = Bu=yd and u = up
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in Feasible branch
%                         >0 = Error in Sufficient branch
%                         -1,1 = Solver error (unbounded solution)
%                         -3,3 = Iteration limit exceeded
%                         -4   = Objective error saturates emax
%         itlim      = Number of iterations remaining after solution found
% 
% Calls:
%        DBcaLP1f_sol  -- Subfunction for feasibility branch
%        DBcaLP1s_sol  -- Subfunction for sufficient  branch
%
% Notes:
%    Feasibility is assessed with a hard-coded tolerance on the cost, J, of 1e-5. This
%        should be reset based on the scaling of the inputs.
%
%    The "preferred" control, up, is used to initialize the feasibility branch. The resulting error
%  components |B*up-yd| are used as slack variables to drive the solution toward yd. An upper
%  limit on the objective error components are needed to pose the problem for the bounded
%  solver, and must necessarily be emax(i) >= abs(B(i,:)*up-yd(i)) for the initial solution
%  to be feasible.Because the simplex reduces cost at each step, a sufficient condition on
% emax is emax(i) >= w'*abs(B*up-yd)/w(i);
%
%    Error code < 0 implies an error in the first branch and there is no guarantee on
%  the quality of the output solution other than the control limits and 
%     wd'*(B*u-yd) <= wd'*(B*yp-yd)
%   
%    Error code > 0 for errors in second branch and, while it may not be the minimal 
%  control error, the resulting solution should have B*u=yd and can be used.
%  
%
% Modification History
%   2002      Roger Beck  Original (DPcaLP1it.m)
%   8/2014    Roger Beck  Update for use in text
%

%Figure out how big problem is (use standard CA definitions for m & n
[n,m] = size(B);

%Call Feasibility branch
[u,J, inBout, eout, errout,itlim] = DBcaLP1f_sol(yd,B,wd,emax,up,uMin,uMax,n,m,itlim);

feas = 0;
%Check if feasible...if so call sufficiency branch if not exit
if J < 1e-5
    feas = 1;
	[u,Js, errout,itlim] = DBcaLP1s_sol(B*u,B,wu,u,up,inBout, eout, uMin,uMax,n,m,itlim);   
    if Js < 1e-5
       feas = 2;
    end
end
return;

end % DB_CALPcaLP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,J, inBout, eout, errout,itlim] = DBcaLP1f_sol(yd,B,w,emax,up,uMin,uMax,n,m,itlim)
%   Dual Branch Control Allocation--Feasibility Branch
%      Objective Error Minimization
%
% [u,J, inBout, eout, errout,itlim] = DBcaLP1f_sol(yd,B,w,emax,up,uMin,uMax,n,m,itlim);
%
%    Determines if a feasible control solution can be found minimizing the error between
%  the desired and obtained objective (i.e. is yd in the AMS?). If yd is unattainable
%  the returned solution minimizes the weighted 1-norm of the objective error.
%  The Bounded Revised Simplex solver is called to minimize
%    min J= | diag(wd)*(B*u - yd) |_1  s.t. umin <= u <=umax
%   (See A.1.2 and Example A.2 in the text for a discussion of a similar formulation).
%
%
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          w [n,1]  = Weighting on Objective error
%          up [m,1]  = Preferred control solution (used for initialization)
%          emax[n,1] = Upper Bound for objective error
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          n         = Number of objectives
%          m         = Number of controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         J          = Cost=weighted 1-norm of objective error wd'*abs(Bu-yd)
%         eout       = Bound flag for controls (true at lower bound, false at upper bound)
%         errout     = Error Status code
%                         0 = found solution
%                        -3 = Iteration limit exceeded in branch
%                        -1 = Solver error
%                        -4 = Objective error for solution has at least one component
%                               at emax.
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%    The "preferred" control, up, is used to initialize the solver and the resulting error
%  components |B*up-yd| are used as slack variables to drive the solution toward yd. An upper
%  limit on the objective error components are needed to pose the problem for the bounded
%  solver, and must necessarily be emax(i) >= abs(B(i,:)*up-yd(i)) for the initial solution
%  to be feasible.Because the simplex reduces cost at each step, a sufficient condition on
% emax is emax(i) >= w'*abs(B*up-yd)/w(i);
%
% Modification History
%   2002      Roger Beck  Original
%   8/2014    Roger Beck  Update
%

%Initialize error code to zero
errout = 0;

%Formulate as an LP problem
A = [eye(n) -eye(n) -B];
b = B*uMin-yd;
c = [w;w;zeros(m,1)];
h = [emax; emax; uMax-uMin];


% A feasible initial condition is the up, using the objective error as
%  slack variables.
eyd = B*up-yd;
x0 = [ max(eyd,zeros(n,1));max(-eyd,zeros(n,1));zeros(m,1)];

%Find Basis Variables for initial solution
%   If preferred control has zero objective error in an axis, identify
%     additional basic variables that are zero.
%   This is unlikely with floating point data and limited precision
%     --could also handle by biasing zero terms in eyp by eps.
%
%
indn = 1:n;
numzer = length(find(eyd==0));
inBi = [indn(eyd>0) n+indn(eyd<0) (2*n):( (2*n)-1+numzer )];

e = true(2*n+m,1);

%Solve using Bounded Revised Simplex
[y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inBi,h,e,n,2*n+m,itlim);

%Construct solution to original LP problem from bounded simplex output
%  Set non-basic variables to 0 or h based on e2
%  Set basic variables to y2 or h-y2.
xout = zeros(2*n+m,1);
xout(inB2) = y2;
xout(~e2) = -xout(~e2)+h(~e2);


if itlim<=0
    errout = -3;
    disp('Too Many Iterations Finding Final Solution');
end
if errsimp
    errout = -1;
    disp('Solver error');
end

%Check if solution contains error terms that are limited at their upper limit
tmp = ~e2;
tmp(inB2) = false;
if any(tmp(1:2*n))
   err = -4;
   disp('Output objective error at bounds');
end
%Compute cost (objective error)-- used to determine feasibility in
%subsequent program
J = c'*xout;

%Convert solution back to control variable
u = xout(2*n+1:2*n+m)+uMin;

%Output controls in basis for subsequent program
%  If solution is feasible, then any error terms in basis are 0
%   and can be replaced with one of the controls at zero

inBout = inB2(inB2>(2*n))-2*n;
if length(inBout) < n  %If there are too few controls in basis
    cind = 1:m;
    cvec = setdiff(cind,inBout);
    inBout = [inBout cvec(1:(n-length(inBout)))];
end

eout = e2(2*n+1:2*n+m);

end %DBcaLP1f_sol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,J,errout,itlim] = DBcaLP1s_sol(yd,B,w,u0,up,inBi, ei, uMin,uMax,n,m,itlim)
%   Dual Branch Control Allocation--Sufficient Branch
%      Control Error Minimization
%
% function [u,J,errout,itlim] = DBcaLP1s_sol(yd,B,w,u0,up,inBi, ei, uMin,uMax,n,m,itlim);
%
%   Assumes that the desired objective, yd, is attainable and seeks to minimize the error
%  between the controls and a "preferred" control solution.    
%  The Bounded Revised Simplex solver is called to minimize
%    min J= | diag(wu)*(u - up) |_1  s.t. umin <= u <=umax and Bu=yd
%   (See A.4.1 in the text for a discussion of a similar formulation).
%
%
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          w [n,1]  = Weighting on Control error
%          u0 [m,1]  = Initial basic solution that attains Bu=yd
%          up [m,1]  = Preferred control solution (used for initialization)
%          inBi[n,1] = List of controls in the basis set(i.e. not at limits)
%          ei[n,1]   = Controls not at upper bounds (ei(i) = false if ui(i) = uMax(i))
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          n         = Number of objectives
%          m         = Number of controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         J          = Cost=weighted 1-norm of control error wd*abs(u-up)
%         errout     = Error Status code
%                         0 = found solution
%                         3 = Iteration limit exceeded in  branch
%                         1 = Solver error
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%    The input solution u0 is assumed to be a basic feasible solution (i.e. yd is attainable
%    and m-n controls are at a limit).
%
% Modification History
%   2002      Roger Beck  Original
%   8/2014    Roger Beck  Update 
%


%Figure out how big problem is (use standard CA definitions for m & n
%[n,m] = size(B);

%Initialize Error Code
errout = 0;

%Formulate as an LP problem
A = [B -B];
b = yd-B*up;
c = [w;w];
h = [uMax-up;up-uMin];

%Any feasible solution that satisfies Bu = yd will satisfy constraint
eup = u0-up;
x0 = [max(eup,zeros(m,1));max(-eup,zeros(m,1))];

%Use identified basis from feasible problem--note that limited variables
%above are not in the basis.

inB = [inBi(eup(inBi)>=0) inBi(eup(inBi)<0)+m];


%Determine which non-basic variables are at upper bound
%
e = (x0 < h/2);
e(inB) = true;

%Solve using Bounded Revised Simplex
%e = true(2*m,1);
%e(evec) = false;

[y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB,h,e,n,2*m,itlim);


%Construct solution to original LP problem from bounded simplex output
%  Set non-basic variables to 0 or h based on e2
%  Set basic variables to y2 or h-y2.
xout = zeros(2*m,1);
xout(inB2) = y2;
xout(~e2) = -xout(~e2)+h(~e2);


if itlim<=0
    errout = 3;
    disp('Too Many Iterations Finding Final Solution');
end
if errsimp
    errout = 1;
    disp('Solver error');
end
%Compute cost output
J = c'*xout;
%Convert solution back to control variable
u = xout(1:m)-xout(m+1:2*m)+up;

end %DBcaLP1s_sol

function [u, feas, errout,itlim] = DBinf_LPCA(yd,B,wd,up,wu,emax,uMin,uMax,itlim) % note
%   Dual Branch Control Allocation - Linear Program
%      Objective Error Minimization (1-norm)
%      Control Error Minimization (inf-norm)
%
% [u, feas, errout,itlim] = DB_LPCA(yd,B,wd,up,wu,emax,uMin,uMax,itlim);
%
%    Uses a Bounded Revised Simplex solver to solve two linear Programming problems
%  The first branch, feasibility, seeks to minimize the weighted 1-norm of the
%    objective error
%    min J= | diag(wd)*(B*u - yd) |_1  s.t. umin <= u <=umax
%
%   If J = 0 then yd is on the interior of the AMS and a second program seeks to minimize
%   the maximum of a set of weighted absolute errors compared to a "preferred" control solution.
%    min  [ max (i=1:n) wu(i)*abs(u(i)-up(i))]  s.t. umin <= u <= umax and  Bu=yd
%
%   (In the text, the overall structure is presented in A.5.1)
%
%   (The overall structure of this allocator is similar to that presented
%   in Buffington, J. et al. "On-Line System Identification for Aircraft with Distributed
%         Control Effectors", Int. J. Robust Nonlinear Control 9, 1033-1049 (1999), however
%    the control error minimization branch is modified).
%        
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          wd [n,1]  = Weighting on Objective error
%          up [m,1]  = Preferred control solution
%          wu [m,1]  = Weighting on control error
%          emax[n,1] = Upper Bound for objective error
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         feas       = Feasible flag
%                        0 = yd is unachievable
%                        1 = Bu=yd
%                        2 = Bu=yd and u = up
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in Feasible branch
%                         >0 = Error in Sufficient branch
%                         -1,1 = Solver error (unbounded solution)
%                         -3,3 = Iteration limit exceeded
%                         -4   = Objective error saturates emax
%         itlim      = Number of iterations remaining after solution found
% 
% Calls:
%        DBcaLP1f_sol  -- Subfunction for feasibility branch
%        DBcaLP1s_sol  -- Subfunction for sufficient  branch
%
% Notes:
%    Feasibility is assessed with a hard-coded tolerance on the cost, J, of 1e-5. This
%        should be reset based on the scaling of the inputs.
%
%    The "preferred" control, up, is used to initialize the feasibility branch. The resulting error
%  components |B*up-yd| are used as slack variables to drive the solution toward yd. An upper
%  limit on the objective error components are needed to pose the problem for the bounded
%  solver, and must necessarily be emax(i) >= abs(B(i,:)*up-yd(i)) for the initial solution
%  to be feasible.Because the simplex reduces cost at each step, a sufficient condition on
% emax is emax(i) >= w'*abs(B*up-yd)/w(i);
%
%    Error code < 0 implies an error in the first branch and there is no guarantee on
%  the quality of the output solution other than the control limits and 
%     wd'*(B*u-yd) <= wd'*(B*yp-yd)
%   
%    Error code > 0 for errors in second branch and, while it may not be the minimal 
%  control error, the resulting solution should have B*u=yd and can be used.
%  
%
% Modification History
%   2002      Roger Beck  Original (DPcaLP1it.m)
%   8/2014    Roger Beck  Update for use in text
%

%Figure out how big problem is (use standard CA definitions for m & n
[n,m] = size(B);

%Call Feasibility branch
[u,J, inBout, eout, errout,itlim] = DBcaLP1f_sol(yd,B,wd,emax,up,uMin,uMax,n,m,itlim);

feas = 0;
%Check if feasible...if so call sufficiency branch if not exit
if J < 1e-5
    feas = 1;
	[u,Js, errout,itlim] = DBinfcaLP1s_sol(yd,B,wu,u,up,inBout, eout, uMin,uMax,n,m,itlim);   
    if Js < 1e-5
       feas = 2;
    end
end
return;

end % DB_CALPcaLP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,J,errout,itlim] = DBinfcaLP1s_sol(yd,B,w,u0,up,inBi, ei, uMin,uMax,n,m,itlim)
%   Dual Branch Control Allocation--Sufficient Branch
%      Control Error Minimization
%
% function [u,J,errout,itlim] = DBcaLP1s_sol(yd,B,w,u0,up,inBi, ei, uMin,uMax,n,m,itlim);
%
%   Assumes that the desired objective, yd, is attainable and seeks to minimize the maximum
%   individual component of the absolute error between the controls and a "preferred" control solution.
%   (See A.4.2 in the text for a discussion of a similar formulation).
%  The Bounded Revised Simplex solver is called to minimize
%    min J= | diag(wu)*(u - up) |_inf  s.t. umin <= u <=umax
%
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          w [n,1]  = Weighting on Control error
%          u0 [m,1]  = Initial basic solution that attains Bu=yd
%          up [m,1]  = Preferred control solution (used for initialization)
%          inBi[n,1] = List of controls in the basis set(i.e. not at limits)
%          ei[n,1]   = Controls not at upper bounds (ei(i) = false if ui(i) = uMax(i))
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          n         = Number of objectives
%          m         = Number of controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         J          = Cost=weighted inf-norm of control error max w(i)*abs(u(i)-up(i))
%         errout     = Error Status code
%                         0 = found solution
%                         3 = Iteration limit exceeded in  branch
%                         1 = Solver error
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%    The input solution u0 is assumed to be a basic feasible solution (i.e. yd is attainable
%    and m-n controls are at a limit).
%
% Modification History
%   2002      Roger Beck  Original
%   8/2014    Roger Beck  Update 
%   8/2015    Roger Beck  Update from 1-norm to inf-norm 
%


%Figure out how big problem is (use standard CA definitions for m & n
%[n,m] = size(B);

%Initialize Error Code
errout = 0;

%Formulate as an LP problem
A = [B -B zeros(n,2*m+1); ...
     diag(w) zeros(m) eye(m) zeros(m) -ones(m,1) ; ...
     zeros(m) diag(w) zeros(m) eye(m) -ones(m,1); ...
     ];
b = [yd-B*up; zeros(2*m,1)];
c = [zeros(4*m,1); 1];
h = [uMax-up;up-uMin; uMax-up;up-uMin; max(max(uMax-up),max(up-uMin))];

%Any feasible solution that satisfies Bu = yd will satisfy constraint
eup = u0-up;
[us,is] = max([eup; -eup]);
x01 = [max(eup,zeros(m,1));max(-eup,zeros(m,1))];
x02 = us-x01;
x0 = [x01;x02;us];

%Use identified basis from feasible problem--note that limited variables
%above are not in the basis.

inBt = 2*m + (1:2*m);
inBt(is) = [];
inB = [inBi(eup(inBi)>=0) inBi(eup(inBi)<0)+m inBt 4*m+1];


%Determine which non-basic variables are at upper bound
%
e = (x0 < h/2);
e(inB) = true;

%Solve using Bounded Revised Simplex
%e = true(2*m,1);
%e(evec) = false;

[y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB,h,e,n,2*m,itlim);


%Construct solution to original LP problem from bounded simplex output
%  Set non-basic variables to 0 or h based on e2
%  Set basic variables to y2 or h-y2.
xout = zeros(4*m+1,1);
xout(inB2) = y2;
xout(~e2) = -xout(~e2)+h(~e2);

if itlim<=0
    errout = 3;
    disp('Too Many Iterations Finding Final Solution');
end
if errsimp
    errout = 1;
    disp('Solver error');
end
%Compute cost output
J = c'*xout;
%Convert solution back to control variable
u = xout(1:m)-xout(m+1:2*m)+up;

end %DBinfcaLP1s_sol

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u, errout] = DP_LPCA(yd,B,uMin,uMax,itlim)
% Direction Preserving Control Allocation Linear Program
%
% function [u, errout] = DP_LPCA(yd,B,uMin,uMax,itlim);
%
%    Solves the control allocation problem while preserving the
%  objective direction for unattainable commands. The solution
%  is found by solving the problem,
%    min -lambda,
%    s.t. B*u = lambda*yd, uMin<=u<=uMax, 0<= lambda <=1
%
%  For yd outside the AMS, the solution returned is that the
%  maximum in the direction of yd.
%
%  For yd strictly inside the AMS, the solution achieves
%  Bu=yd and m-n controls will be at their limits; but there
%  is no explicit preference to which solution will be 
%  returned. This limits the usefulness of this routine as
%  a practical allocator unless preferences for attainable solutions
%  are handled externally.
%
%  (For derivation of a similar formulation see A.1.2 and A.2.3 in the
%  text)
%
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in finding initial basic feasible solution
%                         >0 = Error in finding final solution
%                         -1,1 = Solver error (unbounded solution)
%                         -2   = Initial feasible solution not found
%                         -3,3 = Iteration limit exceeded
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%   If errout ~0 there was a problem in the solution. %
%
%    Error code < 0 implies an error in the initialization and there is no guarantee on
%  the quality of the output solution other than the control limits.
%    Error code > 0 for errors in final solution--B*u is in the correct direction and has
%  magnitude < yd, but B*u may not equal yd (for yd attainable)
%   or be maximized (for yd unattainable)
%
% Modification History
%   2002      Roger Beck  Original (DPcaLP8.m)
%   8/2014    Roger Beck  Update for use in text


%Initialize error code to zero
errout = 0;

tol = 1e-10;
%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

%Check to see if yd == 0
%  May want to adjust the tolerance to improve numerics of later steps
if (all(abs(yd) < tol))    %yd = 0 ==> u=0
    errout = -1;
    u = zeros(m,1);
    return;
end

%Construct an LP using scaling parameter to enforce direction preserving
A = [B -yd];
b = -B*uMin;
c = [zeros(m,1);-1];
h = [uMax-uMin; 1];


%To find Feasible solution construct problem with appended slack variables
sb = 2*(b > 0)-1;
Ai = [A diag(sb)];   
ci = [zeros(m+1,1);ones(n,1)];
inBi = m+2:m+n+1;
ei = true(m+n+1,1);
hi = [h;2*abs(b)];

%Use Bounded Revised Simplex to find initial basic feasible point of
%original program
[y1, inB1, e1,itlim, errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n,m+n+1,itlim);

%Check that Feasible Solution was found
if itlim<=0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1>(m+1))
    errout = -2;
    disp('No Initial Feasible Solution found');
end
	if errsimp
	    errout = -1;
		disp('Solver error');
	end

if errout ~=0  % Construct an incorrect solution to accompany error flags
        xout = zeros(m+1,1);
        indv = inB1<=(m+1);
        xout(inB1(indv)) = y1(indv);
        xout(~e1(1:m+1)) = -xout(~e1(1:m+1))+h(~e1(1:m+1));
else  % No Error continue to solve problem
    
    
    %Solve using initial problem from above
    [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB1,h,e1(1:m+1),n,m+1,itlim);
    
    %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(m+1,1);
    xout(inB2) = y2;
    xout(~e2) = -xout(~e2)+h(~e2);
    
    if itlim<=0
        errout = 3;
        disp('Too Many Iterations Finding Final Solution');
    end
	if errsimp
	    errout = 1;
		disp('Solver error');
	end
    
end


%Transform back to control variables
u = xout(1:m)+uMin;
return;
end

function [u, errout] = DPprio_LPCA(yd,ye,B,uMin,uMax,itlim)
% Direction Preserving Control Allocation Linear Program
%
% function [u, errout] = DP_LPCA(yd,B,uMin,uMax,itlim);
%
%    Solves the control allocation problem while preserving the
%  objective direction for unattainable commands. The solution
%  is found by solving the problem,
%    min -lambda,
%    s.t. B*u = lambda*yd, uMin<=u<=uMax, 0<= lambda <=1
%
%  For yd outside the AMS, the solution returned is that the
%  maximum in the direction of yd.
%
%  For yd strictly inside the AMS, the solution achieves
%  Bu=yd and m-n controls will be at their limits; but there
%  is no explicit preference to which solution will be 
%  returned. This limits the usefulness of this routine as
%  a practical allocator unless preferences for attainable solutions
%  are handled externally.
%
%  (For derivation of a similar formulation see A.1.2 and A.2.3 in the
%  text)
%
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in finding initial basic feasible solution
%                         >0 = Error in finding final solution
%                         -1,1 = Solver error (unbounded solution)
%                         -2   = Initial feasible solution not found
%                         -3,3 = Iteration limit exceeded
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%   If errout ~0 there was a problem in the solution. %
%
%    Error code < 0 implies an error in the initialization and there is no guarantee on
%  the quality of the output solution other than the control limits.
%    Error code > 0 for errors in final solution--B*u is in the correct direction and has
%  magnitude < yd, but B*u may not equal yd (for yd attainable)
%   or be maximized (for yd unattainable)
%
% Modification History
%   2002      Roger Beck  Original (DPcaLP8.m)
%   8/2014    Roger Beck  Update for use in text


%Initialize error code to zero
errout = 0;

tol = 1e-10;
%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

%Check to see if yd == 0
%  May want to adjust the tolerance to improve numerics of later steps
if (all(abs(yd+ye) < tol))    %yd = 0 ==> u=0
    errout = -1;
    u = zeros(m,1);
    return;
end

%Construct an LP using scaling parameter to enforce direction preserving
A = [B -yd];
b = ye-B*uMin;
c = [zeros(m,1);-1];
h = [uMax-uMin; 1];


%To find Feasible solution construct problem with appended slack variables
sb = 2*(b > 0)-1;
Ai = [A diag(sb)];   
ci = [zeros(m+1,1);ones(n,1)];
inBi = m+2:m+n+1;
ei = true(m+n+1,1);
hi = [h;2*abs(b)];

%Use Bounded Revised Simplex to find initial basic feasible point of
%original program
[y1, inB1, e1,itlim, errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n,m+n+1,itlim);

%Check that Feasible Solution was found
if itlim<=0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1>(m+1))
    errout = -2;
    disp('No Initial Feasible Solution found');
end
	if errsimp
	    errout = -1;
		disp('Solver error');
	end

if errout ~=0  % Construct an incorrect solution to accompany error flags
    if all( abs(ye)<=tol )
        xout = zeros(m+1,1);
        indv = inB1<=(m+1);
        xout(inB1(indv)) = y1(indv);
        xout(~e1(1:m+1)) = -xout(~e1(1:m+1))+h(~e1(1:m+1));
    else
        [u,~] = DPprio_LPCA(ye,[0;0;0],B,uMin,uMax,itlim);
%         [u,~] = SBprio_LPCA(ye,[0;0;0],B,w,up,uMin,uMax,itlim);
        return;
    end
else  % No Error continue to solve problem
    
    
    %Solve using initial problem from above
    [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB1,h,e1(1:m+1),n,m+1,itlim);
    
    %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(m+1,1);
    xout(inB2) = y2;
    xout(~e2) = -xout(~e2)+h(~e2);
    
    if itlim<=0
        errout = 3;
        disp('Too Many Iterations Finding Final Solution');
    end
	if errsimp
	    errout = 1;
		disp('Solver error');
	end
    
end


%Transform back to control variables
u = xout(1:m)+uMin;
return;
end

function [u,itlim,errout] = DPscaled_LPCA(yd,B,uMin,uMax,itlim)
% Direction Preserving Control Allocation Linear Program
%     Reduced formulation (Solution Scaled from Boundary)
%
% function [u,itlim,errout] = DPscaled_LPCA(yd,B,uMin,uMax,itlim);
%
%    Solves the control allocation problem while preserving the
%  objective direction for unattainable commands. The reduced
%  dimension of the linear program passed to the Bounded Revised
%  Simplex solver is formed by forcing the solution to be on the
%  boundary of the AMS and eliminating the highest magnitude
%  objective by solving the other constraints in terms of it.
%
%  For yd outside the AMS, the solution returned is that the
%  maximum in the direction of yd
%    B*u= lamda*yd
%    max lamda s.t. uMin <= u <= uMax
%
%  Reducing the degrees of freedom elminates the problems of redundant
%  solutions for attainable objectives. If the desired objective is on the
%  interior of the AMS the solution is scaled from the solution on the
%  boundary, yielding the same controls as the Direct Allocation solution.
%  
%  (In the text this solution is discussed in section A.5.3)
%
%   (See Bodson, M., "Evaluation of Optimization Methods for
%          Control Allocation",  AIAA 2001-4223).
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in finding initial basic feasible solution
%                         >0 = Error in finding final solution
%                         -1,1 = Solver error (unbounded solution)
%                         -2   = Initial feasible solution not found
%                         -3,3 = Iteration limit exceeded
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%    If yd is close to zero, u = 0;
%
%    Error code < 0 implies an error in the initialization and there is no guarantee on
%  the quality of the output solution other than the control limits.
%    Error code > 0 for errors in final solution.
%
% Modification History
%   2002      Roger Beck  Original ( DPcaLP2.m)
%   8/2014    Roger Beck  Update


%Initialize error code to zero
errout = 0;

tol = 1e-10;
%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

% Locate the maximum magnitude element in the desired objective
[my,iy]=max(abs(yd));

%Trivial solution, if desired moment is close to zero
%  May want to adjust the tolerance to improve numerics of later steps
if (my < tol)    %yd = 0 ==> u=0
    errout = -1;  %Set flag to let caller know that it wasn't solved
    u = zeros(m,1);
    return;
end

%Transform Problem by Reordering Objectives with maximum first
Bt = B([iy setdiff(1:n,iy)],:);
ydt = yd([iy setdiff(1:n,iy)]);
ydt(2:3) = ydt([3 2]);Bt([2 3],:) = Bt([3 2],:);
%%Convert into a LP problem
M = [ydt(2:n) -ydt(1)*eye(n-1)];
A = M*Bt;
b = -A*uMin;
c = -Bt'*ydt;
h = uMax-uMin;

%To find Feasible solution construct problem with appended slack variables
sb = 2*(b > 0)-1;
Ai = [A diag(sb)];
ci = [zeros(m,1);ones(n-1,1)];
inBi = m+1:m+n-1;
ei = true(m+n-1,1);
hi = [h;2*abs(b)];
%Use Bounded Revised Simplex to find initial basic feasible point
[y1, inB1, e1,itlim,errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n-1,m+n-1,itlim);

%Check that Feasible Solution was found
if itlim<=0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1>m)
    errout = -2;
    disp('No Initial Feasible Solution found');
end
if errsimp
    errout = -1;
    disp('Solver error');
end

if errout ~=0  % Construct an incorrect solution to accompany error flags
%     xout = zeros(m,1);
%     xout(inB1(1:m)) = y1(1:m);
%     xout(~e1(1:m)) = -xout(~e1(1:m))+h(~e1(1:m));
    xout = zeros(m,1);
    indv = inB1<=(m);
    xout(inB1(indv)) = y1(indv);
    xout(~e1(1:m)) = -xout(~e1(1:m))+h(~e1(1:m));
    
else  % No Error continue to solve problem
    
    
    
    %Solve using initial problem from above
    [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB1,h,e1(1:m),n-1,m,itlim);
    
    %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(m,1);
    xout(inB2) = y2;
    xout(~e2) = -xout(~e2)+h(~e2);
    
    if itlim<=0
        errout = 3;
        disp('Too Many Iterations Finding Final Solution');
    end
    if errsimp
        errout = 1;
        disp('Solver error');
    end
    
    
end

%Transform Solution Back Into control variables
% u(i) = x(i)+umin(i) if e(i)
u = xout+uMin;
%Rescale controls so solution is not on boundary of Omega.
rho = ydt'*Bt*u/(ydt'*ydt);
if rho > 1
    u = u/rho;
end


return;
end

function [u,errout] = MO_LPCA(yd,B,up,lam, eMax,uMin,uMax,itlim)
% Mixed Optimization (Single Branch) Control Allocation Linear Program
%    Objective Error Minimizing
%    Control Error minimizing
%
% function [u,errout] = MO_LPCA(yd,B,w,up,uMin,uMax,itlim);
%
%    Solves the control allocation problem while seeking to
%  simultaneously minimize the error in the desired objective
%  and the error in the controls using a single linear program.
%
%  Finds the solution that minimizes
%  min |B*u - yd |_1+lambda*|(u-up)|_1
%   such that
%   uMin <= u <= uMax
%
%  The balance between the two objectives is determined by the
%  weight lambda. This weight needs to be chosen small so that
%  the objective error minimization part of the problem dominates
%  the solution in order to ensure that B*u is close to yd for
%  attainable solutions.
%
%  (Section A.5.2 and example A.7 in the text discuss Single Branch
%   optimization routines including this formulation).
%
%   (See Bodson, M., "Evaluation of Optimization Methods for
%          Control Allocation",  AIAA 2001-4223).
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          up[m,1]   = Preferred Control Vector
%          lam       = Control Error Weighting (single parameter)
%          eMax[n,1] = Maximum Objective error
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                         0 = found solution
%                         -1 = Solver error (unbounded solution)
%                         -3 = Iteration limit exceeded
%                         -4 = Initial feasible solution saturates emax
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%   Attainable controls may result in a small error in objective depending on 
%   scaling of control and objective errors and the magnitude of the weight, lambda.
%
%   Solution with errors, errout~=0, only ensures that control limits are respected and
%  that weighted error, |B*u - yd |_1+lambda*|(u-up)|_1, is <= |B*up-yd|_1.
%
%   The "preferred" control, up, is used to initialize the optimization. The resulting error
%  components |B*up-yd| are used as slack variables to drive the solution toward yd. An upper
%  limit on the objective error components are needed to pose the problem for the bounded
%  solver, and must necessarily be emax(i) >= abs(B(i,:)*up-yd(i)) for the initial solution
%  to be feasible. Because the simplex reduces cost at each step, a sufficient condition on
% emax is emax(i) >= w'*abs(B*up-yd)/w(i);
%
%   If errout ~0 there was a problem in the solution. %
%
% Modification History
%   2002      Roger Beck  Original (MOcaLP1.m(
%   8/2014    Roger Beck  Update for use in text

%Initialize error code to zero
errout = 0;

%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

%Formulate as an LP problem
A = [eye(n) -eye(n) -B B];
b = B*up-yd;
c = [ones(2*n,1); lam*ones(2*m,1)];
h = [eMax;eMax;uMax-up;up-uMin];

% A feasible initial condition is the up, using the objective error as
%  slack variables.
eyp = B*up-yd;
x0 = [ max(eyp,zeros(n,1));max(-eyp,zeros(n,1));zeros(2*m,1)];

%Find Basis Variables for initial solution
%   If preferred control has zero objective error in an axis, identify
%     additional basic variables that are zero.
%   This is unlikely with floating point data and limited precision
%     --could also handle by biasing zero terms in eyp by eps.
%  
%
indn = 1:n;
numzer = length(find(eyp==0));
inBi = [indn(eyp>0) n+indn(eyp<0) (2*n):( (2*n)-1+numzer )];
e = true(2*(n+m),1);

%Solve using Bounded Revised Simplex
    [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inBi,h,e,n,2*(m+n),itlim);

   %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(2*(m+n),1);
    xout(inB2) = y2;
    xout(~e2) = -xout(~e2)+h(~e2);
       
    %Check if solution contains error terms that are limited at their upper limit
tmp = ~e2;
tmp(inB2) = false;
if any(tmp(1:2*n))
   errout = -4;
   disp('Output objective error at bounds');
end

    if itlim<=0
        errout = -3;
        disp('Too Many Iterations Finding Final Solution');
    end
    if errsimp
        errout = -1;
        disp('Solver error');
    end

%%Transform Solution Back Into control variables
% Note that x(2*n+1:2*n+m) are the + differences from the preferred control
% solution and x(2*n+m+1:2*(n+m)) are the - differences.

u = xout(2*n+1:2*n+m)-xout(2*n+m+1:2*(n+m))+up;
return;
end

function [u,errout] = SB_LPCA(yd,B,w,up,uMin,uMax,itlim) % note
% Single Branch Control Allocation Linear Program
%    Direction Preserving
%    Control Error minimizing
%
% function [u,errout] = SB_LPCA(yd,B,w,up,uMin,uMax,itlim);
%
%    Solves the control allocation problem while seeking to
%  simultaneously preserve the direction for unattainable
%  objectives and minimizing the control error for attainable
%  commands using a single linear program.
%
%  Finds the solution that minimizes
%  min -lambda + |diag(w)*(u-up)|_1
%   such that
%  B*u = lambda*yd
%   uMin <= u <= uMax
%      0 <= lambda <= 1
%
%  The balance between the two objectives is determined by the
%  weight vector, w. For unattainable moments the constraints
%  ensure that the direction of the command is maintained, 
%  however, the weights should be small compared to the
%  relative scaling of the units on the control vector and the
%  objective vector, so that the error term doesn't pull the
%  optimum solution away from B*u =yd.
%
%  (Section A.5.2 and example A.6 in the text discuss Single Branch
%  optimization routines of including this formulation).
% 
%   (See Buffington, J. "Tailess Aircraft Control Allocation",
%      AIAA-97-3695 for a similar approach that seeks to minimize
%      control usage (i.e. up = 0) and partitions lambda to prioritize
%      command components)
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          w [m,1]   = Control Error Weighting
%          up[m,1]   = Preferred Control Vector
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in finding initial basic feasible solution
%                         >0 = Error in finding final solution
%                         -1,1 = Solver error (unbounded solution)
%                         -2   = Initial feasible solution not found
%                         -3,3 = Iteration limit exceeded
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%   Attainable controls may result in a small error in objective depending on 
%   scaling of control and objective errors and the magnitude of the weights.
%
%    Error code < 0 implies an error in the initialization and there is no guarantee on
%  the quality of the output solution other than the control limits.
%    Error code > 0 for errors in final solution, result has B*u in the right direction
% and magnitude <= yd.
%
% Modification History
%   2002      Roger Beck  Original (SBcaLP2)
%   8/2014    Roger Beck  Update for use in text


%Initialize error code to zero
errout = 0;

%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

%Formulate as an LP problem
A = [B -B -yd];
b = -B*up;
c = [w; w; -1];
h = [uMax-up; up-uMin;1];

%To find Feasible solution construct problem with appended slack variables
sb = 2*(b > 0)-1;
Ai = [A diag(sb)];   
ci = [zeros(2*m+1,1);ones(n,1)];
inBi = 2*m+2:2*m+n+1;
ei = true(2*m+n+1,1);
hi = [h;2*abs(b)];

%Use Bounded Revised Simplex to find initial basic feasible point
[y1, inB1, e1,itlim,errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n,2*m+n+1,itlim);

%Check that Feasible Solution was found
if itlim<=0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1>(2*m+1))
    errout = -2;
    disp('No Initial Feasible Solution found');
end
if errsimp
    errout = -1;
    disp('Solver error');
end

if errout ~=0  % Construct an incorrect solution to accompany error flags
    xout = zeros(2*m+1,1);
    indv = inB1<=(2*m+1);
    xout(inB1(indv)) = y1(indv);
    xout(~e1(1:2*m+1)) = -xout(~e1(1:2*m+1))+h(~e1(1:2*m+1));
    
    
else  % No Error continue to solve problem
    
        
    %Solve using initial problem from above
    [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB1,h,e1(1:(2*m+1)),n,(2*m+1),itlim);
    
    %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(2*m+1,1);
    xout(inB2) = y2;
    xout(~e2) = -xout(~e2)+h(~e2);
    
    if itlim<=0
        errout = 3;
        disp('Too Many Iterations Finding Final Solution');
    end
    if errsimp
        errout = 1;
        disp('Solver error');
    end

    
end
    
%Transform Solution Back Into control variables
% Note that x(1:m) are the + differences from the preferred control
% solution and x(m+1:m) are the - differences.

u = xout(1:m)-xout(m+1:2*m)+up;
return;
end

function [u,errout] = SBprio_LPCA(yd,ye,B,w,up,uMin,uMax,itlim) % note
% Single Branch Control Allocation Linear Program
%    Direction Preserving
%    Control Error minimizing
%
% function [u,errout] = SB_LPCA(yd,B,w,up,uMin,uMax,itlim);
%
%    Solves the control allocation problem while seeking to
%  simultaneously preserve the direction for unattainable
%  objectives and minimizing the control error for attainable
%  commands using a single linear program.
%
%  Finds the solution that minimizes
%  min -lambda + |diag(w)*(u-up)|_1
%   such that
%  B*u = lambda*yd
%   uMin <= u <= uMax
%      0 <= lambda <= 1
%
%  The balance between the two objectives is determined by the
%  weight vector, w. For unattainable moments the constraints
%  ensure that the direction of the command is maintained, 
%  however, the weights should be small compared to the
%  relative scaling of the units on the control vector and the
%  objective vector, so that the error term doesn't pull the
%  optimum solution away from B*u =yd.
%
%  (Section A.5.2 and example A.6 in the text discuss Single Branch
%  optimization routines of including this formulation).
% 
%   (See Buffington, J. "Tailess Aircraft Control Allocation",
%      AIAA-97-3695 for a similar approach that seeks to minimize
%      control usage (i.e. up = 0) and partitions lambda to prioritize
%      command components)
%
%  Inputs:
%          yd [n]    = Desired objective
%          B [n,m]   = Control Effectiveness matrix
%          w [m,1]   = Control Error Weighting
%          up[m,1]   = Preferred Control Vector
%          uMin[m,1] = Lower bound for controls
%          uMax[m,1] = Upper bound for controls
%          itlim     = Number of allowed iterations limit
%                         (Sum of iterations in both branches)
%
% Outputs:
%         u[m,1]     = Control Solution
%         errout     = Error Status code
%                         0 = found solution
%                         <0 = Error in finding initial basic feasible solution
%                         >0 = Error in finding final solution
%                         -1,1 = Solver error (unbounded solution)
%                         -2   = Initial feasible solution not found
%                         -3,3 = Iteration limit exceeded
%         itlim      = Number of iterations remaining after solution found
%
% Calls:
%         simplxuprevsol = Bounded Revised Simplex solver (simplxuprevsol.m)
%
% Notes:
%   Attainable controls may result in a small error in objective depending on 
%   scaling of control and objective errors and the magnitude of the weights.
%
%    Error code < 0 implies an error in the initialization and there is no guarantee on
%  the quality of the output solution other than the control limits.
%    Error code > 0 for errors in final solution, result has B*u in the right direction
% and magnitude <= yd.
%
% Modification History
%   2002      Roger Beck  Original (SBcaLP2)
%   8/2014    Roger Beck  Update for use in text

%Tolerance for unknown == 0
tol = 1e-10;

%Initialize error code to zero
errout = 0;

%Figure out how big the problem is (use standard CA definitions for m & n)
[n,m] = size(B);

%Formulate as an LP problem
A = [B -B -yd];
b = -B*up+ye;
c = [w; w; -1];
h = [uMax-up; up-uMin;1];

%To find Feasible solution construct problem with appended slack variables
sb = 2*(b > 0)-1;
Ai = [A diag(sb)];   
ci = [zeros(2*m+1,1);ones(n,1)];
% inBi = [2*m+2:2*m+n+1];
inBi=false(1,2*m+n+1);
inBi(2*m+2:2*m+n+1) =1;
ei = true(2*m+n+1,1);
hi = [h;2*abs(b)];

%Use Bounded Revised Simplex to find initial basic feasible point
% [y1, inB1, e1,itlim,errsimp] = simplxuprevsol(Ai,ci',b,inBi,hi,ei,n,2*m+n+1,itlim);
[y1, inB1, e1,itlim,errsimp] = simpl(Ai,ci',b,inBi,hi,ei,n,2*m+n+1,itlim);
%Check that Feasible Solution was found
if itlim<=0
    errout = -3;
    disp('Too Many Iterations Finding initial Solution');
end
if any(inB1>(2*m+1))
% if any(find(inB1)>(2*m+1))
    errout = -2;
    disp('No Initial Feasible Solution found');
end
if errsimp
    errout = -1;
    disp('Solver error');
end

if errout ~=0  % Construct an incorrect solution to accompany error flags
    if all( abs(ye)<=tol )
        xout = zeros(2*m+1,1);
        indv = inB1<=(2*m+1);
        xout(inB1(indv)) = y1(indv);
        xout(~e1(1:2*m+1)) = -xout(~e1(1:2*m+1))+h(~e1(1:2*m+1));
    else
%         [u,~,~] = DPscaled_LPCA(ye,B,uMin,uMax,itlim);
        [u,~] = SBprio_LPCA(ye,[0;0;0],B,w,up,uMin,uMax,itlim);
        return;
    end
else  % No Error continue to solve problem
    
        
    %Solve using initial problem from above
%     [y2, inB2, e2,itlim,errsimp] = simplxuprevsol(A ,c',b,inB1,h,e1(1:(2*m+1)),n,(2*m+1),itlim);
    inBi=false(1,2*m+1);
    inBi(inB1) =1;
    [y2, inB2, e2,itlim,errsimp] = simpl(A,c',b,inBi,h,e1(1:(2*m+1)),n,2*m+1,itlim);
    %Construct solution to original LP problem from bounded simplex output
    %  Set non-basic variables to 0 or h based on e2
    %  Set basic variables to y2 or h-y2.
    xout = zeros(2*m+1,1);
    xout(inB2) = y2;
    xout(~e2) = -xout(~e2)+h(~e2);
    
    if itlim<=0
        errout = 3;
        disp('Too Many Iterations Finding Final Solution');
    end
    if errsimp
	    errout = 1;
		disp('Solver error');       
    end

    
end
    
%Transform Solution Back Into control variables
% Note that x(1:m) are the + differences from the preferred control
% solution and x(m+1:m) are the - differences.

u = xout(1:m)-xout(m+1:2*m)+up;
return;
end

function [y0, inB, e,itlim,errout] = simpl(A,ct,b,inBx,h,e,varargin)
%  Bounded Revised Simplex
%
%function [yout, inBout,bout, itout,errout] = simplxuprevsol(A,ct,b,inB,inD,h,e,m,n,itlim)
%
%   Solves the linear program:
%          minimize c'y 
%          subject to 
%          Ay = b
%          0<= y <= h
%
%  Inputs: 
%          A [m,n]   = lhs Matrix of equaltity constraints
%          ct [1,n]  = transpose of cost vector
%          b [m,1]   = rhs vector for equality constraint
%-----------------------------------------------------------------------
%          inB [m]   = Vector of indices of unknowns in the initial basic set
%          inD [n-m] = Vector of indices of unknowns not in the initial basic set
% from inBx [n], inBx is a logical vector that non-zeros element indices
% are inB, so as inD
%--------------------------------------------------------------------------------
%          h[n,1]    = Upper Bound for unknowns
%          e[n,1]    = Sign for unknown variables (+ lower bound, - upper bound)
%  Optional inputs:
%          m,n       = number of constraints, unknowns (Opposite standard
%                      CA convention
%          itlim     = Upper bound on the allowed iterations\
%
% Outputs:
%         yout[n,1]  = Optimal output variable
%         inBout     = indices of Basic vectors in output
%         eout       = sign associate with output unknowns
%         itout      = number of iterations remaining out of itlim
%         errout     = Flag (=true) if unbounded is set
%
% Modification History
%   2002      Roger Beck  Original
%   8/2014    Roger Beck  Update for use
%   9/2014    Roger Beck  Added anti-cycling rule

%Optional Inputs
switch  length(varargin)
    case 0
    itlim = inf;
    [m,n] = size(A);
    case 1
    itlim = varargin{1};
    [m,n] = size(A);
    case 2
    itlim = inf;
    m = varargin{1};
    n = varargin{2};
    case 3
    itlim =varargin{3};
    m = varargin{1};
    n = varargin{2};
 end    	
    	
%Tolerance for unknown == 0
tol = 1e-10;

%Index list for non-basic variables
nind = 1:(n-m);

%Partition A
inB=find(inBx);
inD=find(~inBx);
%-----------------------
% inD = setdiff(1:n, inB);
%-----------------------

%Adjust signs problem if variables are initialized at upper
% bounds.
A(:,~e) = -A(:,~e);
ct(~e) = -ct(~e);
b = b + A(:,~e)*h(~e);

y0 = A(:,inB)\b;  %Initial Solution

%Initialize Loop Termination Conditions
done = false;
unbounded = false;

%Main Simplex loop
while (~done  || ~unbounded ) && (itlim > 0)
    itlim = itlim-1;

    %Calculate transpose of relative cost vector based on current basis
    lamt = ct(inB)/A(:,inB);
    rdt = ct(inD)-lamt*A(:,inD);
    %Find minimum relative cost
    [minr, qind] = min(rdt);
    if minr >=0  % If all relative costs are positive then the solution is optimal
        done = true;
        break;
    end
    qel = inD(qind);  % Unknown to Enter the basis minimizes relative cost
    yq = A(:,inB)\A(:,qel); %Vector to enter in terms of the current Basis vector
    
    if all(abs(yq)<=tol)
      unbounded = true;
      disp(' Solution is unbounded');  % Check this condition
      break
    end

    %Compute ratio how much each current basic variable will have to move for the entering
    % variable.

    rat = y0./yq; 
    
    % If yq < 0 then increasing variable when it leaves the basis will minimize cost
    hinB = h(inB);
    indm = yq<0;
    rat(indm) = rat(indm) - hinB(indm)./yq(indm);
    % If an element yq ~=0 then it doesn't change for the entering variable and shouldn't
    %  be chosen
    indz = abs(yq)<=tol;
    rat(indz) = inf;

    % Variable to exit is moving to its minimum value
    [minrat, p] = min(rat);

   % If the minimum ratio is zero, then the solution is degenerate and the entering
   %   variable will not change the basis---invoke Bland's selection rule to avoid
   %   cycling.
    if (abs(minrat) <= tol)
       % Find negative relative cost
       indm = nind(rdt<0); %Note that since minr <0 indm is not empty   
       qind = indm(1);
       qel = inD(qind);  % Unknown to Enter the basis is first indexed to avoid cycling
       yq = A(:,inB)\A(:,qel); %Vector to enter in terms of the current Basis vector
       if all(abs(yq)<=tol)
           unbounded = true;
           disp(' Solution is unbounded');  % Check this condition
           break
       end
       % Recompute rations and determine variable to leave
       rat = y0./yq; 
        % If yq < 0 then increasing variable when it leaves the basis will minimize cost
        hinB = h(inB);
        indm = yq<0;
        rat(indm) = rat(indm) - hinB(indm)./yq(indm);
        % If an element yq ~=0 then it doesn't change for the entering variable and shouldn't
        %  be chosen
        indz = abs(yq)<=tol;
        rat(indz) = inf;

        % Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
        [minrat, p] = min(rat);
    end

  % Maintain the bounded simplex as only having lower bounds by recasting 
  % any variable that needs to move to its opposite bound.
    if (minrat >= h(qel))
           %Case 1: Entering variable goes to opposite bound and current basis is maintained
            e(qel) = ~e(qel);
            A(:,qel) = -A(:,qel);
             b = b + A(:,qel)*h(qel);
             ct(qel) = -ct(qel);
    elseif yq(p) > 0
           %Case 2: Leaving variable returns to lower bound (0)	
           pel = inB(p);
           inB(p)= qel;
           inD(qind)= pel;
     else
           %Case 2: Leaving variable moves to upper bound	
            pel = inB(p);
            e(pel)=~e(pel);
            A(:,pel) = -A(:,pel);
            inB(p)= qel;
            inD(qind)= pel;
            ct(pel) = -ct(pel);
            b = b + A(:,pel)*h(pel);
     end
        
    y0 = A(:,inB)\b; % Compute new Basic solution;
end
errout = unbounded;     
end

function [y0, inB, e,itlim,errout] = simplxuprevsol(A,ct,b,inB,h,e,varargin)
%  Bounded Revised Simplex
%
%function [yout, inBout,bout, itout,errout] = simplxuprevsol(A,ct,b,inB,inD,h,e,m,n,itlim)
%
%   Solves the linear program:
%          minimize c'y 
%          subject to 
%          Ay = b
%          0<= y <= h
%
%  Inputs: 
%          A [m,n]   = lhs Matrix of equaltity constraints
%          ct [1,n]  = transpose of cost vector
%          b [m,1]   = rhs vector for equality constraint
%          inB [m]   = Vector of indices of unknowns in the initial basic set
%          inD [n-m] = Vector of indices of unknowns not in the initial basic set
%          h[n,1]    = Upper Bound for unknowns
%          e[n,1]    = Sign for unknown variables (+ lower bound, - upper bound)
%  Optional inputs:
%          m,n       = number of constraints, unknowns (Opposite standard
%                      CA convention
%          itlim     = Upper bound on the allowed iterations\
%
% Outputs:
%         yout[n,1]  = Optimal output variable
%         inBout     = indices of Basic vectors in output
%         eout       = sign associate with output unknowns
%         itout      = number of iterations remaining out of itlim
%         errout     = Flag (=true) if unbounded is set
%
% Modification History
%   2002      Roger Beck  Original
%   8/2014    Roger Beck  Update for use
%   9/2014    Roger Beck  Added anti-cycling rule

%Optional Inputs
switch  length(varargin)
    case 0
    itlim = inf;
    [m,n] = size(A);
    case 1
    itlim = varargin{1};
    [m,n] = size(A);
    case 2
    itlim = inf;
    m = varargin{1};
    n = varargin{2};
    case 3
    itlim =varargin{3};
    m = varargin{1};
    n = varargin{2};
 end    	
    	
%Tolerance for unknown == 0
tol = 1e-8;

%Index list for non-basic variables
nind = 1:(n-m);

%Partition A
inD = setdiff(1:n, inB);

%Adjust signs problem if variables are initialized at upper
% bounds.
A(:,~e) = -A(:,~e);
ct(~e) = -ct(~e);
b = b + A(:,~e)*h(~e);

y0 = A(:,inB)\b;  %Initial Solution

%Initialize Loop Termination Conditions
done = false;
unbounded = false;

%Main Simplex loop
while (~done  || ~unbounded ) && (itlim > 0)
    itlim = itlim-1;

    %Calculate transpose of relative cost vector based on current basis
    lamt = ct(inB)/A(:,inB);
    rdt = ct(inD)-lamt*A(:,inD);
    %Find minimum relative cost
    [minr, qind] = min(rdt);
    if minr >=0  % If all relative costs are positive then the solution is optimal
        done = true;
        break;
    end
    qel = inD(qind);  % Unknown to Enter the basis minimizes relative cost
    yq = A(:,inB)\A(:,qel); %Vector to enter in terms of the current Basis vector
    
    if all(abs(yq)<=tol)
      unbounded = true;
      disp(' Solution is unbounded');  % Check this condition
      break
    end

    %Compute ratio how much each current basic variable will have to move for the entering
    % variable.

    rat = y0./yq; 
    
    % If yq < 0 then increasing variable when it leaves the basis will minimize cost
    hinB = h(inB);
    indm = yq<0;
    rat(indm) = rat(indm) - hinB(indm)./yq(indm);
    % If an element yq ~=0 then it doesn't change for the entering variable and shouldn't
    %  be chosen
    indz = abs(yq)<=tol;
    rat(indz) = inf;

    % Variable to exit is moving to its minimum value
    [minrat, p] = min(rat);

   % If the minimum ratio is zero, then the solution is degenerate and the entering
   %   variable will not change the basis---invoke Bland's selection rule to avoid
   %   cycling.
    if (abs(minrat) <= tol)
       % Find negative relative cost
       indm = nind(rdt<0); %Note that since minr <0 indm is not empty   
       qind = indm(1);
       qel = inD(qind);  % Unknown to Enter the basis is first indexed to avoid cycling
       yq = A(:,inB)\A(:,qel); %Vector to enter in terms of the current Basis vector
       if all(abs(yq)<=tol)
           unbounded = true;
           disp(' Solution is unbounded');  % Check this condition
           break
       end
       % Recompute rations and determine variable to leave
       rat = y0./yq; 
        % If yq < 0 then increasing variable when it leaves the basis will minimize cost
        hinB = h(inB);
        indm = yq<0;
        rat(indm) = rat(indm) - hinB(indm)./yq(indm);
        % If an element yq ~=0 then it doesn't change for the entering variable and shouldn't
        %  be chosen
        indz = abs(yq)<=tol;
        rat(indz) = inf;

        % Variable to exit is moving to its minimum value--Note that min returns the lowest index minimum
        [minrat, p] = min(rat);
    end

  % Maintain the bounded simplex as only having lower bounds by recasting 
  % any variable that needs to move to its opposite bound.
    if (minrat >= h(qel))
           %Case 1: Entering variable goes to opposite bound and current basis is maintained
            e(qel) = ~e(qel);
            A(:,qel) = -A(:,qel);
             b = b + A(:,qel)*h(qel);
             ct(qel) = -ct(qel);
    elseif yq(p) > 0
           %Case 2: Leaving variable returns to lower bound (0)	
           pel = inB(p);
           inB(p)= qel;
           inD(qind)= pel;
     else
           %Case 2: Leaving variable moves to upper bound	
            pel = inB(p);
            e(pel)=~e(pel);
            A(:,pel) = -A(:,pel);
            inB(p)= qel;
            inD(qind)= pel;
            ct(pel) = -ct(pel);
            b = b + A(:,pel)*h(pel);
     end
        
    y0 = A(:,inB)\b; % Compute new Basic solution;
end
errout = unbounded;     
end
