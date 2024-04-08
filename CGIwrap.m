function [u] = CGIwrap(IN_MAT,v,NumU)
% Make single input, single output version of Cascading Generalized Inverse
% control allocation for use in Simulink via the MATLAB Fcn block
% IN_MAT = [B     d
%           umin' 0
%           umax' 0
%           INDX  0]
% 20140524  KAB Modified to include INDX, which is used to specify active
%               effectors
% 20160316  KAB Modified to simplify implementation
% global NumU
% Get sizes
[k2,m1]=size(IN_MAT);
k=k2-3;
m=m1-1;
% If matrices too small, set contols to zero and return
if k<1 || m<1 || norm(IN_MAT)<1e-16
    u=zeros(NumU,1);
    return
end
% Partition input matrix into component matrices
B=IN_MAT(1:k,1:m);
% v=IN_MAT(1:k,end);
umin=IN_MAT(k+1,1:m)';
umax=IN_MAT(k+2,1:m)';
INDX=IN_MAT(k+3,1:m)';

% get active effectors
B_act=B(:,INDX>0.5); 
umin_act=umin(INDX>0.5);
umax_act=umax(INDX>0.5);

[n,m]=size(B_act);
[u_act]=CGIfunc(v,B_act,umin_act,umax_act);
u=zeros(NumU,1);
u(INDX>0.5,1)=u_act;

end


function [u] = CGIfunc(d_in,Bin,umin_in,umax_in)
% Single output version

% Define variables for later use
% various "small" values used to check if approximately zero
eps_col=0.001; % small number to check if B column is too small
eps_lim=0.001; % small number to check if limits collapsed
eps_row=0.001; % small number to check if B row is too small
eps_ans=0.001; % small number to prevent divide by zero when computing scale
n=0; % number of rows = number of controlled vairables
m=0; % number of columns = number of active effectors
B=[];
col_ok=0;
[nin,m_in]=size(Bin);
hitlim=0;
done=0;
sing=0;
pmax=0;
d0=d_in;
forcep=0;

% Check columns of Bin and control limits
% Build B from valid columns of Bin
for j=1:m_in
    % Check that limits have not collapsed
    if(umax_in(j)-umin_in(j))>eps_lim
        col_ok=0;
        % Check that column of B > (small number)
        if(max(abs(Bin(:,j))) > eps_col)
            col_ok=1;
        end %if
        % If it's good, use it
        if(col_ok == 1)
            m=m+1;
            indx(m,1)=j;
            umin(m,1)=umin_in(j);
            umax(m,1)=umax_in(j);
            B(:,m)=Bin(:,j);
            u(j,1)=0.0;
        else
            % If it's no good, set u as average of limits
            u(j,1)=0.5*(umax_in(j)+umin_in(j));
        end %if
    else
        % If limits collapsed, set u as average of limits
        u(j,1)=0.5*(umax_in(j)+umin_in(j));
    end %if
end % for j

% Account for effects of controls removed by changing d
d_in=d_in-Bin*u;

% This is needed because both limits may be +'ve or -'ve
% This assumes umax>umin, hopefully this will always be true
for i=1:m
    if(umax(i)<0&&umin(i)<0)
        ut2(i,1)=umax(i);
    elseif(umin(i)>0&&umax(i)>0)
        ut2(i,1)=umin(i);
    else
        ut2(i,1)=0;
    end % if
end % for j
if(m>nin)
    B1=B;
    m1=m;
    indx1=indx;
    umin1=umin;
    umax1=umax;
else
    m1=m;
end

passes=0;

%% From here on down, repeat until solution
% is found or no controls left

while done==0 && m>=1
    passes=passes+1;
    Pt2=zeros(m_in,nin);
    % Remove rows of B if necessary
    BB=[];
    d=[];
    n=0;
    for i=1:nin
        row_ok=0;
        if(max(abs(B(i,:)))>eps_row)
            row_ok=1;
        end %if
        % If row of B > (small number) use it
        if(row_ok==1)
            n=n+1;
            d(n,1)=d_in(i);
            BB(n,:)=B(i,:);
        end % if
    end % for i
    % If no rows left, done
    if n==0;
        done=1;
    else
        B=BB;
        % Center the problem
        %   Try to redo centering to preserve direction
        uc=[];
        uc=ut2;
        uminc=umin-uc;
        umaxc=umax-uc;
        dc=d-B*uc;
        
        % Compute gain scaling matrix based on initial off-set between d_in
        % and dc
        if passes==1;
            dscale=eye(n,n);
            for i=1:n
                if abs(d0(i)) > 0.1
                    dscale(i,i)=dc(i)/d0(i);
                end % if
            end %for i
        end %if
        % Get controls, solve B*ut=dc for ut
        if forcep==1;
            dc=zeros(size(dc));
            forcep=0;
        end %if
        
        [ut,Pt,infot]=get_u3(dc,B);
            singt=infot(1);
            pmaxt=infot(2);
            pindx=infot(3);
            
            if singt > sing
                sing=singt;
            end %if
            
            if forcep==0
                % if pmaxt > pmax
                pmax = pmaxt;
                % end; %if
            else
                forcep=0;
            end %if
            
            % Check limits
            hitlim=0;
            scale=1.0;
            imax=0;
            for i=1:m
                temp=1.0;
                if(ut(i)<uminc(i))
                    hitlim=1;
                    temp=uminc(i)/min(ut(i),-eps_ans);
                elseif(ut(i)>umaxc(i))
                    hitlim=1;
                    temp=umaxc(i)/max(ut(i),eps_ans);
                end; % if
                if(temp<scale)
                    scale=temp;
                    imax=i;
                end %if
            end %for i
            
            % Remove controls as necessary
            if(hitlim==1)
                % Set max limited control to limiting value
                % and recompute d_in to accout for its effects
                u(indx(imax))=uc(imax)+(ut(imax)*scale);
                d_in=d_in-Bin(:,indx(imax))*u(indx(imax));
                % Redo indx
                % handle scaling properly
                us=uc+ut.*scale;
                ut2=[];
                ut2=[us(1:imax-1);us(imax+1:m)];
                indx2=[indx(1:imax-1);indx(imax+1:m)];
                m=m-1;
                indx=indx2;
                
                % Redo umin,umax,up,Wp,B,ut2
                if(m~=0)
                    umin=[];
                    umax=[];
                    B=[];
                    
                    for i=1:m
                        umin(i,1)=umin_in(indx(i));
                        umax(i,1)=umax_in(indx(i));
                        B(:,i)=Bin(:,indx(i));
                    end %for i
                end %if
            else
                % If none hit limits, done
                % Undo the centering stuff
                ut2=[];
                ut2=ut+uc;
                done=1;
            end %if
    end %if
end %while
% If done, reassemble control vector
for i=1:m
    u(indx(i))=ut2(i);
end %for i

info(1)=sing;
info(2)=passes;
info(3)=pmax;

end %function CGIfunc

% function [u,P,info]=get_u3(d,B,Wpi,Wd,PL,PL_n)
function [u,P,info]=get_u3(d,B)
% Solves Bu=d
% Uses unweighted Moore-Penrose pseudo-inverse 
% d=nx1
% B=nxm
% u=mx1
%
% program checks that all matrices are correct size
% if not, u=zeros(up)
% sing=1 means matrix was singular

% Check Matrix Sizes
[dn,d1]=size(d);
[n,m]=size(B);

if dn~=n||d1~=1
    u=zeros(m,1);
    info(1)=1;
    info(2)=1;
    P=zeros(m,n);
    a='Error in get_u:Matrices Wrong Size!'
    return
end %if

% 3 cases to check

if m > n
    % Use Minimum Norm with Weights and Prefs
    % Minimizes weighted norm of control vector (weighted |u|)
    %
    dum=(B*(B'));
    [dumi,sing]=get_inv(dum);
    P=B'*dumi;
elseif m==n
    % Use inverse on square matrix
    [P,sing]=get_inv(B);
else
    % Use Least-Squares
    % Minimize error |Bu-d|
    %
    btb=B'*B;
    [btbi,sing]=get_inv(btb);
    P=btbi*B';
end %if

if sing==0
    % Matrix that was inverted was not singular
    u=P*d;
elseif m>=2
    a='singular matrix!'
    % Combine into single effector
    BS=B(:,1);
    for i=2:m
        BS=BS+B(:,i);
    end %for i
    % use least-squares
    % Minimize error |Bu-d|
    %
    btb=(BS')*BS;
    [btbi,sing2]=get_inv(btb);
    P=btbi*(BS');
    if sing2==0
    % Matrix that was inverted was not singular
        u=P*d;
    else
    % Matrix that was inverted was singular
    % Use zeros for effectors and P matrix
        u=zeros(m,1);
        P=zeros(m,n);
    end %if
else
    u=zeros(m,1);
    P=zeros(m,n);
end %if

info(1)=sing;
info(2)=0;
info(3)=0;

end %function get_u3

function [AI,sing]=get_inv(A)
% [AI,sing]=get_inv(A)
% Inverts A using simple algebra
% sing = 1 if matrix is singular
%
% Warning!! Matrix must be square and 3x3 or smaller!

eps=0.0001; %if abs(determinant) greater than this, treat as non-singular
sing=0;

[n,m]=size(A);

if n~=m
    AI=eye(n,m);
    sing=1;
    a='Error1 in get_inv:Matrix Not Square!'
elseif n==3
    adj11=A(2,2)*A(3,3)-A(2,3)*A(3,2);
    adj12=-1*(A(2,1)*A(3,3)-A(2,3)*A(3,1));
    adj13=A(2,1)*A(3,2)-A(2,2)*A(3,1);
    adj21=-1*(A(1,2)*A(3,3)-A(1,3)*A(3,2));
    adj22=A(1,1)*A(3,3)-A(1,3)*A(3,1);
    adj23=-1*(A(1,1)*A(3,2)-A(1,2)*A(3,1));
    adj31=A(1,2)*A(2,3)-A(1,3)*A(2,2);
    adj32=-1*(A(1,1)*A(2,3)-A(1,3)*A(2,1));
    adj33=A(1,1)*A(2,2)-A(1,2)*A(2,1);
    det=adj11*A(1,1)+adj12*A(1,2)+adj13*A(1,3);
    if abs(det) > eps
        deti=1/det;
        AI(1,1)=adj11*deti;
        AI(1,2)=adj21*deti;
        AI(1,3)=adj31*deti;
        AI(2,1)=adj12*deti;
        AI(2,2)=adj22*deti;
        AI(2,3)=adj32*deti;
        AI(3,1)=adj13*deti;
        AI(3,2)=adj23*deti;
        AI(3,3)=adj33*deti;
    else
        AI=zeros(size(A));
        sing=1;
    end %if
elseif n==2
    det=A(1,1)*A(2,2)-A(1,2)*A(2,1);
    if abs(det) > eps
        deti=1/det;
        AI(1,1)=A(2,2)*deti;
        AI(1,2)=-A(1,2)*deti;
        AI(2,1)=-A(2,1)*deti;
        AI(2,2)=A(1,1)*deti;
    else
        AI=zeros(size(A));
        sing=1;
    end %if
elseif n==1
    det=A(1,1);
    if abs(det) > eps
        AI(1,1)=1/det;
    else
        AI(1,1)=0;
        sing=1;
    end %if
else
    AI=eye(n,m);
    sing=1;
    a='Error2 in get_inv:Matrix Wrong Size!'
end %if
end %function get_inv
