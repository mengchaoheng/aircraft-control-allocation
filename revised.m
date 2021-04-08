classdef revised < handle
    % REVISED This class implements a Revised Simplex Method to solve a
    % linear programming problem in the following format
    %   min/max c'x
    %   s.t.   Ax {>=, =, <=} b
    %   x >= 0
    %
    % The code here follows the variation taught in the FuO Course for
    % OR Msc. Specifically, this class is designed for class demonstration
    % and small problems only, neither suitable for large problem nor
    % powerful enough for high performance computing.
    %
    % Example: See example.m
    %
    % 8 October 2013
    % Yiming Yan
    % University of Edinburgh
    
    %% Properties
    % Problem data
    properties (SetAccess = private)
        A;                  % Porb data: coefficient matrix for constraints
        b;                  % Prob data: rhs
        c;                  % Prob data: coefficient for objective function
        
        inq;                % m dimentional vector,
                            % 1 for >=; 0 for ==; -1 for <=.
        
        m;                  % Number of constraints
        n;                  % Number of variables
        fval;               % Optimal objective value
        type = 'min';       % Min or max type
    end
    
    % Iteration information
    properties (SetAccess = private)
        x;					% Vector of primal variables
        y;					% Vector of dual variables
        z;                  % Current objective value
        vars_list;          % List of variables' names
        
        d;					% Vector of reduced costs
        
        dNq;                % The value of the most negative reduced cost
        q; 					% Nonbasic variable to enter the basis
        t;					% Basic variable to leave the basis
        
        
        updCol;				% Updated column
        
        sigma;				% Step size
        basis;              % Current basis
        nonbasis;			% Current nonbasis
        M;                  % Big M
        counter;            % Iteration counter
        terminate;          % flag for termination: 1 terminate; 0 continue
        status;             % maxIter, unbounded, infeasible, optimal
    end
    
    % Constants
    properties( Constant )
        maxIter = 20;       %  max number of iterations allowed
    end
    
    
    %% Methods
    methods
        %% Constructor
        function p = revised(A,b,c,inq,minmax)
            % Constructor Read in given problem data including:
            %   (A,b,c) : Problem data
            %      inq  : A vector specifiees the inequalties and equaliteis
            %             1 for >=; 0 for ==; -1 for <=.
            %   minmax  : specifies 'min' or 'max'
            
            % Check input
            if nargin < 4
                error('REVISED: Not enough input');
            end
            
            %% Read the input
            % Check min or max
            if strcmpi(minmax,'max')
                c = -c;
                p.type = 'max';
            end
            
            % Check if b is negative
            b_Neg      = b < 0;
            inq(b_Neg) = -inq(b_Neg);
            b(b_Neg)   = -b(b_Neg);
            A(b_Neg,:) = -A(b_Neg,:);
            
            % Assiqn properties
            p.A = A; p.b = b; p.c = c;
            p.inq = inq;
            p.terminate = 0; p.counter = 0;
            
            % Set the value of big M
            p.M = ceil(max(norm(A),norm(b))/100)*100;
            
            [p.m, p.n] = size(A);
            
            p.basis = zeros(1,p.m);
            % Form the vars_list
            p.vars_list = cell(1,p.n);
            for i=1:p.n
                p.vars_list{i} = ['x' num2str(i)];
            end
        end
        
        %% Driver function
        function solve(p)
            % SOLVE This function actually solves the given problem
            % using a revised simplex method.
            
            p.transform_to_standardForm;
            p.initialise;
            % Main loop
            while 1
                p.increment_counter;
                p.compute_dual_vars;                  p.printInfo('BTRAN');       
                p.price_nonbasic_vars;                p.printInfo('PRICE');       
                p.choose_nonbasis_var_to_enter_basis; p.printInfo('ChuzC');    
                p.find_updated_col;                   p.printInfo('FTRAN');
                p.find_basic_var_to_leave_basis;      p.printInfo('ChuzR');
                p.update_step;
                p.update_basis;                       p.printInfo('Update');
                if p.terminate
                    break;
                end
            end
            p.output_summary;
        end
    end
    
    %% Internal functions
    methods (Access = private)
        function transform_to_standardForm(p)
            % add slack variables
            n_slack = sum(p.inq < 0) + sum(p.inq > 0);
            p.A = [p.A zeros(p.m, n_slack)];
            p.c = [p.c; zeros(n_slack,1)];
            
            idx_slack = p.n + 1;
            for i = 1:p.m
                switch p.inq(i)
                    case 1  % >=
                        p.A(i, idx_slack) = -1;
                        idx_slack = idx_slack + 1;
                        p.vars_list{end+1} = ['s' num2str(i)];
                    case -1 % <=
                        p.A(i, idx_slack) = 1;
                        idx_slack = idx_slack + 1;
                        p.vars_list{end + 1} = ['s' num2str(i)];
                end
            end
            
            p.update_dimension;
        end
        
        
        function initialise(p)
            % Find potential basis
            for i = 1:p.n
                if nnz(p.A(:,i)) == 1
                    row_number = find(p.A(:,i) == 1);
                    if ~isempty(row_number)
                        % Add the col if the current row has not been selected
                        if p.basis(row_number) == 0
                            p.basis(row_number) = i;
                        end
                    end
                end
            end
            
            % Add bigM if necessary
            n_artfVar = p.m - sum(p.basis > 0);
            if  n_artfVar > 0
                % Record the index for artificial variables
                p.A = [p.A zeros(p.m, n_artfVar)];
                
                add_to_rows = find(p.basis == 0);
                
                % Formulate matrix A with newly added artificial vars
                for i = 1:length(add_to_rows)
                    p.A(add_to_rows(i), p.n + i) = 1;
                    p.vars_list{end+1} = ['a' num2str(i)];
                    p.basis(add_to_rows(i)) = p.n + i;
                end
                
                p.c = [p.c; p.M*ones(n_artfVar,1)];
            end
            p.nonbasis = setdiff(1:p.m,p.basis);
            p.update_dimension;
            
            p.x = zeros(p.n,1);
            p.x(p.basis) =  p.A(:,p.basis)\p.b;
            
            p.d = p.c - p.A'*p.c(p.basis);
            p.z = 0;
        end
        
        function compute_dual_vars(p)
            p.y = inv( p.A(:,p.basis) )'*p.c(p.basis);
        end
        
        function price_nonbasic_vars(p)
            p.d(p.nonbasis) = p.c(p.nonbasis) - p.A(:,p.nonbasis)'*p.y;
            p.d(p.basis) = 0;
        end
        
        function choose_nonbasis_var_to_enter_basis(p)
            p.dNq = min(p.d(p.nonbasis));
            if p.dNq >= 0
                p.terminate = 1;
                p.status = 'optimal';
                if strcmpi(p.type, 'max')
                    p.fval = - p.z;
                else
                    p.fval = p.z;
                end
            else
                p.q = find( p.d(p.nonbasis) == p.dNq );
                p.q = p.nonbasis(p.q);
            end
            
        end
        
        function find_updated_col(p)
            p.updCol = p.A(:,p.basis)\p.A(:,p.q);
        end
        
        function find_basic_var_to_leave_basis(p)
            [p.sigma, indx] = min(p.x(p.basis)./p.updCol);
            if isinf(p.sigma)
                p.terminate = 1;
                p.status = 'unbounded';
            end
            
            p.t = p.basis(indx);
        end
        
        function update_step(p)
            e = zeros(p.n,1); e(p.q) = 1;
            p.x(p.basis) = p.x(p.basis) - p.sigma*p.updCol;
            p.x(p.nonbasis) = p.x(p.nonbasis) + p.sigma*e(p.nonbasis);
            p.z = p.z + p.dNq*p.sigma;
        end
        
        function update_basis(p)
            p.basis(p.basis == p.t) = p.q;
            p.nonbasis(p.nonbasis == p.q) = p.t;
        end
        
        function increment_counter(p)
            p.counter = p.counter + 1;
            if p.counter == p.maxIter
                p.terminate = 1;
                p.status = 'maxIter';
            end
        end
        
        function update_dimension(p)
            p.n = size(p.A,2);
        end
        
       
        function printInfo(p, procedName)
           % PrintInfo Print out iterative information
           % procedName: BTRAN, PRICE, CHUZC, FTRAN, CHUZR, UPDATE
           if ~p.terminate
               digits = 3;
               switch lower(procedName)
                   case 'btran'
                       fprintf('\n=========== Iteration %s ===========\n\n',...
                           num2str(p.counter));
                       fprintf('  B  = %s;\n', mat2str(p.A(:,p.basis)));
                       fprintf('  N  = %s;\n', mat2str(p.A(:,p.nonbasis)));
                       fprintf('  cB = %s;\n', mat2str(p.c(p.basis)));
                       fprintf('  cN = %s;\n\n', mat2str(p.c(p.nonbasis)));
                       

                       fprintf('* BTRAN : y^{T}     = c_{B}^{T}B^{-1}    = %s\n',...
                           mat2str(p.y',digits));
                       
                   case 'price'
                       fprintf('* PRICE : d_{N}^{T} = c_{N}^{T} - y^{T}N = %s\n',...
                           mat2str(p.d(p.nonbasis)',digits));
                       
                   case 'chuzc'
                       fprintf('* ChuzC : Choose the most negative reduced cost (%s) and increase %s.\n',...
                           num2str(p.dNq, digits), p.vars_list{p.q});
                       
                   case 'ftran'
                       fprintf('* FTRAN : Find the column of %s in the updated tableau, and the RHS\n\n',...
                           p.vars_list{p.q});
                       fprintf('\t B^{-1}N_{%s}     = %s\n',...
                           p.vars_list{p.q}, mat2str(p.updCol, digits));
                       fprintf('\t RHS (= B^{-1}b)  = %s\n\n',...
                           mat2str(p.x(p.basis), digits));
                       
                   case 'chuzr'
                       fprintf('* ChuzR : Find the max value of %s that maintains feasiblity\n\n',...
                           p.vars_list{p.q});
                       for i = 1:p.m+2
                           if i==1
                               fprintf('\t Basis |%5s |%5s ||%10s\n',...
                                   p.vars_list{p.q}, '=',  'Limit')
                               fprintf('\t ----------------------------------\n');
                           elseif i == p.m+2
                               fprintf('\t ----------------------------------\n');
                               fprintf('\t       |%5s |%5s ||%10s\n',...
                                   num2str(p.dNq, digits),' ', num2str(p.z));
                           else
                               for j =1:4
                                   tmp = p.x(p.basis);
                                   if j == 1
                                       fprintf('\t %5s |',...
                                           p.vars_list{p.basis(i-1)});
                                   elseif j == 2
                                       tmp_str = num2str(p.updCol(i-1), digits);
                                       if p.basis(i-1) == p.t
                                           tmp_str = [tmp_str '*'];
                                       end
                                       fprintf('%5s |', tmp_str)
                                   elseif j == 3
                                       fprintf('%5s ||',...
                                           num2str(tmp(i-1), digits));
                                   elseif j == 4
                                       fprintf('%10s \n',...
                                           [num2str(tmp(i-1), digits) ' / '...
                                           num2str(p.updCol(i-1), digits)] );
                                   end
                               end % end for j=1:4
                           end  % end if i == 1 
                       end  % end for 1 = 1:p.m + 2
                       
                   case 'update'
                       fprintf('* Update: Increase %s by %s, then\n\n',...
                           p.vars_list{p.q}, num2str(p.sigma,digits));
                       fprintf('\tx = %s, z = %s.\n',...
                           mat2str(p.x, digits), num2str(p.z,digits));
                       fprintf('        %s enters the basis and %s leaves.\n',...
                           p.vars_list{p.q}, p.vars_list{p.t} );
               end
           end
        end
        
        function output_summary(p)
            if p.terminate
                fprintf('\n========== ========== ==========\n');
                fprintf('Terminated with status [ %s ]. \n', p.status);
                switch lower( p.status )
                    case 'optimal'
                        fprintf('Optimal solution: \n');
                        for i=1:p.n
                            fprintf('%s = %5.3f\n', p.vars_list{i}, p.x(i));
                        end
                        fprintf('Optimal objective function value: %5.2f\n', p.fval);
                        
                    case 'unbounded'
                        fprintf('This problem is unbounded. \n');
                        fprintf('Unbounded variables: ');
                        
                    case 'infeasible'
                        fprintf('This problem is infeasible.\n');
                    case 'maxiter'
                        fprintf('Maximum number of iterations reached.\n');
                end
            end
        end
    end
end
