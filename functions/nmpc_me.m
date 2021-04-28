function [t, x, u] = nmpc_me(runningcosts, ...
              constraints, ...
              linearconstraints, system, ...
              mpciterations, N, T, tmeasure, xmeasure, u0, x_ref,xp_measure, ...
              Non_linear_system, param, discret_system_matrices, ...
              printHeader, printClosedloopData, plotTrajectories, ...
              Pedest_prediction, Pedest_dynamics)

% Computes the closed loop solution for the NMPC problem defined by
% the functions
%   runningcosts:         evaluates the running costs for state and control
%                         at one sampling instant.
%                         The function returns the running costs for one
%                         sampling instant.
%          Usage: [cost] = runningcosts(t, x, u)
%                 with time t, state x and control u
%   terminalcosts:        evaluates the terminal costs for state at the end
%                         of the open loop horizon.
%                         The function returns value of the terminal costs.
%          Usage: cost = terminalcosts(t, x)
%                 with time t and state x
%   constraints:          computes the value of the restrictions for a
%                         sampling instance provided the data t, x and u
%                         given by the optimization method.
%                         The function returns the value of the
%                         restrictions for a sampling instance separated
%                         for inequality restrictions c and equality
%                         restrictions ceq.
%          Usage: [c,ceq] = constraints(t, x, u)
%                 with time t, state x and control u
%   terminalconstraints:  computes the value of the terminal restrictions
%                         provided the data t, x and u given by the
%                         optimization method.
%                         The function returns the value of the
%                         terminal restriction for inequality restrictions
%                         c and equality restrictions ceq.
%          Usage: [c,ceq] = terminalconstraints(t, x)
%                 with time t and state x
%   linearconstraints:    sets the linear constraints of the discretized
%                         optimal control problem. This is particularly
%                         useful to set control and state bounds.
%                         The function returns the required matrices for
%                         the linear inequality and equality constraints A
%                         and Aeq, the corresponding right hand sides b and
%                         beq as well as the lower and upper bound of the
%                         control.
%          Usage: [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
%                 with time t, state x and control u
%   system:               evaluates the difference equation describing the
%                         process given time t, state vector x and control
%                         u.
%                         The function returns the state vector x at the
%                         next time instant.
%          Usage: [y] = system(t, x, u, T)
%                 with time t, state x, control u and sampling interval T
% for a given number of NMPC iteration steps (mpciterations). For
% the open loop problem, the horizon is defined by the number of
% time instances N and the sampling time T. Note that the dynamic
% can also be the solution of a differential equation. Moreover, the
% initial time tmeasure, the state measurement xmeasure and a guess of
% the optimal control u0 are required.
%
% Arguments:
%   mpciterations:  Number of MPC iterations to be performed
%   N:              Length of optimization horizon
%   T:              Sampling interval
%   tmeasure:       Time measurement of initial value
%   xmeasure:       State measurement of initial value
%   u0:             Initial guess of open loop control
% Internal Functions:
%   measureInitialValue:          measures the new initial values for t0
%                                 and x0 by adopting values computed by
%                                 method applyControl.
%                                 The function returns new initial state
%                                 vector x0 at sampling instant t0.
%   applyControl:                 applies the first control element of u to
%                                 the simulated process for one sampling
%                                 interval T.
%                                 The function returns closed loop state
%                                 vector xapplied at sampling instant
%                                 tapplied.
%   shiftHorizon:                 applies the shift method to the open loop
%                                 control in order to ease the restart.
%                                 The function returns a new initial guess
%                                 u0 of the control.
%   solveOptimalControlProblem:   solves the optimal control problem of the
%                                 horizon N with sampling length T for the
%                                 given initial values t0 and x0 and the
%                                 initial guess u0 using the specified
%                                 algorithm.
%                                 The function returns the computed optimal
%                                 control u, the corresponding value of the
%                                 cost function V as well as possible exit
%                                 flags and additional output of the
%                                 optimization method.
%   costfunction:                 evaluates the cost function of the
%                                 optimal control problem over the horizon
%                                 N with sampling time T for the current
%                                 data of the optimization method t0, x0
%                                 and u.
%                                 The function return the computed cost
%                                 function value.
%   nonlinearconstraints:         computes the value of the restrictions
%                                 for all sampling instances provided the
%                                 data t0, x0 and u given by the
%                                 optimization method.
%                                 The function returns the value of the
%                                 restrictions for all sampling instances
%                                 separated for inequality restrictions c
%                                 and equality restrictions ceq.
%   computeOpenloopSolution:      computes the open loop solution over the
%                                 horizon N with sampling time T for the
%                                 initial values t0 and x0 as well as the
%                                 control u.
%                                 The function returns the complete open
%                                 loop solution over the requested horizon.
%   dynamic:                      evaluates the dynamic of the system for
%                                 given initial values t0 and x0 over the
%                                 interval [t0, tf] using the control u.
%                                 The function returns the state vector x
%                                 at time instant tf as well as an output
%                                 of all intermediate evaluated time
%                                 instances.
%   printSolution:                prints out information on the current MPC
%                                 step, in particular state and control
%                                 information as well as required computing
%                                 times and exitflags/outputs of the used
%                                 optimization method. The flow of
%                                 information can be controlled by the
%                                 variable iprint and the functions
%                                 printHeader, printClosedloopData and
%                                 plotTrajectories.


    tol_opt = 1e-6;
    % options for fmincon
    options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 10000,...
                'Algorithm', 'active-set',...
                'FinDiffType', 'forward',...
                'RelLineSrchBnd', [],...
                'RelLineSrchBndDuration', 1,...
                'TolConSQP', 1e-6);
    iprint = 5;
    t = []; % Will be used to plot the results
    x = [];
    u = [];
    xp = []; % position of pedestrian
    u_last = [0;0]; % Applied input of the previous time step
    
    
    %%%%
    delta_t = param.dt;
    steps = param.steps;
    %%%%
    
    % Start of the NMPC iteration
    mpciter = 0;
    while(mpciter < mpciterations)
        % Step (1) of the NMPC algorithm:
        %   Obtain new initial value
        t0 = tmeasure;
        x0 = xmeasure;
        
        %%%
        % Compute all matrices of the dynamics for different timescales
        A_d = cell(steps);       
        B_d = cell(steps);
        %fc = cell(steps); % No need for this as f is same for all
        %timescales
        
        for i=1:steps
            % update timescale
            T = delta_t(i);
            % update Ai and Bi
            [A_i, B_i] = discret_system_matrices(t0, x0, [0;0], T, param);
            % Add the new matrices to a new cell
            A_d{i} = A_i;
            B_d{i} = B_i;
        end
         linear.A = A_d;
         linear.B = B_d;
        %%%%%%
          f = Non_linear_system(x0, [0;0], param);  
          linear.f = f;
        %%% Dynamics for First timescale used to
                %update the current state after prediction)
                % compute the constraints 
        linearisation.A = A_d{1};
        linearisation.B = B_d{1};
        linearisation.f = f;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Predict Pedesterian path
        Lane = param.Lane;
        xp0 = xp_measure;
        Xp = Pedest_prediction(T,xp0,N);
    
        % Check if Pedestrian is near the crossing
        s_break = 0; % 
        if (any(Xp(2,:)>= (Lane(1)-param.safety)))&&(any(Xp(2,:)<= (Lane(2)+param.safety)))
            s_break = param.s_break;
        end
    
        % Update Pedestrian position   
        xp = [ xp; xp_measure ];
        [xp_measure,~] = Pedest_dynamics(xp0, T);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 
        % Step (2) of the NMPC algorithm:
        %   Solve the optimal control problem
        t_Start = tic;
        [u_new, V_current, exitflag, output] = solveOptimalControlProblem ...
            (runningcosts, constraints, ...
             linearconstraints, system, ...
            N, t0, x0, u0, T,tol_opt, options, x_ref(:,mpciter+1:mpciter+N+1), ...
        param , linear, u_last,linearisation,s_break);
        t_Elapsed = toc( t_Start );
        
        % compute the predicted trajectory of the car
        x_pred = dynamics(system,N,x0,t0,u_new, param, linear);
        
        %   Print solution : I will check it later
        
        if ( iprint >= 1 )
            printSolution(printHeader, printClosedloopData, ...
                          plotTrajectories, mpciter, x0, u_new, ...
                          iprint, ...
                          exitflag, output, t_Elapsed, param, x_ref,x_pred, xp0, Xp);

        end


        %   Store closed loop data
        t = [ t; tmeasure ];
        x = [ x; xmeasure ]; % CHANGE: I am not using those so I might delete them
        u = [ u; u_new(:,1)];
        %   Prepare restart
        u0 = shiftHorizon(u_new);
        
        % Step (3) of the NMPC algorithm:
        %   Apply control to process
        [tmeasure, xmeasure] = applyControl(system, T, t0, x0, u_new, param, linearisation);
        mpciter = mpciter+1

       u_last = u_new(:,1); 
    end
end

function [tapplied, xapplied] = applyControl(system, T, t0, x0, u, param, linearisation)
    xapplied = system(t0, x0, u(:,1), T, x0, param, linearisation);
    tapplied = t0+T;
end


function [u, V, exitflag, output] = solveOptimalControlProblem ...
    (runningcosts, constraints, ...
    linearconstraints, system, N, t0, x0, u0, T, tol_opt, ...
    options, x_ref, param, linear, u_last,linearisation,x_break)
    x = zeros(length(x0), N+1);
  %  x = computeOpenloopSolution(system, N, T, t0, x0, u0, param, linearisation);
    x = dynamics(system,N,x0,t0,u0, param, linear);
    % Set control and linear bounds
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    for k=1:N
        [Anew, bnew, Aeqnew, beqnew, lbnew, ubnew] = ...
               linearconstraints(t0+k*T,x(:,k),u0(:,k));
        A = blkdiag(A,Anew);
        b = [b, bnew];
        Aeq = blkdiag(Aeq,Aeqnew);
        beq = [beq, beqnew];
        lb = [lb, lbnew];
        ub = [ub, ubnew];
    end

    % Solve optimization problem
    [u, V, exitflag, output] = fmincon(@(u) costfunction(runningcosts, ...
        system, N, T, t0, x0, ...
        u,x_ref,param, linear, u_last), u0, A, b, Aeq, beq, lb, ...
        ub, @(u) nonlinearconstraints(constraints, ...
        system, N, T, t0, x0, u, param, linearisation,x_ref,x_break,linear), options);
end

function x = computeOpenloopSolution(system, N, T, t0, x0, u, param, linearisation)
    x(:,1) = x0;
    for k=1:N
        x(:,k+1) = system(t0, x(:,k), u(:,k), T, x0, param, linearisation);                       
    end
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))];
end

function cost = costfunction(runningcosts, system, ...
                    N, T, t0, x0, u, x_ref, param, ...
                    linear, u_last)
    cost = 0;
    n = length(x0);
    x = zeros(n, N+1);
    % x = computeOpenloopSolution(system, N, T, t0, x0, u, param, linearisation); 
    x = dynamics(system,N,x0,t0,u, param, linear); 
    cost = cost+ runningcosts(t0+1*T, x(:,2), [u_last,u(:,1)], x_ref(:,2), param); % K=0 in Matlab is k=1 and so u(-1) in Matlab u(0)=u_last
    for k=2:(N-1)
        cost = cost+runningcosts(t0+k*T, x(:,k+1), u(:,k-1:k), x_ref(:,k+1), param); % changed to include u(k-1) and x_ref
    end
end

function [c,ceq] = nonlinearconstraints(constraints, ...
    system, ...
    N, T, t0, x0, u, param, linearisation,x_ref, x_break,linear)
    x = zeros(N+1, length(x0));
    %x = computeOpenloopSolution(system, N, T, t0, x0, u, param, linearisation);
    %%%%%%%%
    x = dynamics(system,N,x0,t0,u, param, linear);
    %%%%%%%
    c = [];
    ceq = [];
    for k=1:N
        [cnew, ceqnew] = constraints(t0+k*T,x(:,k),x_ref, x_break, param);
        c = [c cnew];
        ceq = [ceq ceqnew];
    end
    % terminal constraints
    [cnew, ceqnew] = constraints(t0+(N+1)*T,x(:,N+1),x_ref, x_break, param);
    c = [c cnew];
    ceq = [ceq ceqnew];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% OUTPUT %%%%%%%%%%%%%
function printSolution(printHeader, printClosedloopData, ...
             plotTrajectories, mpciter, x0, u, iprint, exitflag, output, t_Elapsed, ...
             param, x_ref, x, xp0, Xp)
    if (mpciter == 0)
        printHeader();
    end
    printClosedloopData(mpciter, u, x0, t_Elapsed);
    switch exitflag
        case -2
        if ( iprint >= 1 && iprint < 10 )
            fprintf(' Error F\n');
        elseif ( iprint >= 10 )
            fprintf(' Error: No feasible point was found\n')
        end
        case -1
        if ( iprint >= 1 && iprint < 10 )
            fprintf(' Error OT\n');
        elseif ( iprint >= 10 )
            fprintf([' Error: The output function terminated the',...
                     ' algorithm\n'])
        end
        case 0
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning IT\n');
        elseif ( iprint >= 10 )
            fprintf([' Warning: Number of iterations exceeded',...
                     ' options.MaxIter or number of function',...
                     ' evaluations exceeded options.FunEvals\n'])
        end
        case 1
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' \n');
        elseif ( iprint >= 10 )
            fprintf([' First-order optimality measure was less',...
                     ' than options.TolFun, and maximum constraint',...
                     ' violation was less than options.TolCon\n'])
        end
        case 2
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning TX\n');
        elseif ( iprint >= 10 )
            fprintf(' Warning: Change in x was less than options.TolX\n')
        end
        case 3
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning TJ\n');
        elseif ( iprint >= 10 )
            fprintf([' Warning: Change in the objective function',...
                     ' value was less than options.TolFun\n'])
        end
        case 4
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning S\n');
        elseif ( iprint >= 10 )
            fprintf([' Warning: Magnitude of the search direction',...
                     ' was less than 2*options.TolX and constraint',...
                     ' violation was less than options.TolCon\n'])
        end
        case 5
        if ( iprint == 1 )
            fprintf('\n');
        elseif ( iprint >= 2 && iprint < 10 )
            fprintf(' Warning D\n');
        elseif ( iprint >= 10 )
            fprintf([' Warning: Magnitude of directional derivative',...
                     ' in search direction was less than',...
                     ' 2*options.TolFun and maximum constraint',...
                     ' violation was less than options.TolCon\n'])
        end
    end
    if ( iprint >= 5 )
    % NOTE: it simply plots the next state from the current x0 wich will be updated with each control point mpciter
        plotTrajectories(x, x0, param, x_ref, xp0, Xp) 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  NON UNIFORM NPMC  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = dynamics(system, N, x0, t0, u0, param, linear)
%% Calculate the state vector X 
% Compute the state vector x by Concatenating all state vectors
% corresponding to different timescales

    %parameters
    A_d = linear.A;    % Matrix A for the 1st timescale (bicycle model)
    B_d = linear.B;    % Matrix B for the 1st timescale
    delta_t = param.dt;
    steps = param.steps;
    linearisation.f = linear.f;
    X = x0;
    Np = N/steps; %length of prediction timescale
    for i=1:steps
        % update timescale
        T = delta_t(i);
        % update Ai and Bi
        linearisation.A = A_d{i};
        linearisation.B = B_d{i};
    
        % Calculate timescale corresponding state vecor
        xi = computeOpenloopSolution(system, Np, T, t0, x0, u0, param, linearisation);
        % update parameters
        x0 = xi(:,Np+1); 
        xi = xi(:,2:Np+1); % remove x(:,1)=x0
        %Np = Np*(i+1); % this is wrong as Np is the same for all timescales
        % Concatenate new vectors
    X = [X,xi];
     
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function [x, t_intermediate, x_intermediate] = dynamic(system, T, t0, x0, u, x_curr)
% % t_intermediate, x_intermediate will be needed to plot the results
%         x = system(t0, x0, u, T, x_curr);
%         x_intermediate = [x0; x];
%         t_intermediate = [t0, t0+T];
% 
% end