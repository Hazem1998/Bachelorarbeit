function [t, x, u] = nmpc_me(runningcosts, ...
              constraints, ...
              linearconstraints, system, ...
              mpciterations, N, T, xmeasure, u0, x_ref,xp_measure,iprint, ...
              Non_linear_system, param, discret_system_matrices, ...
              printHeader, printClosedloopData, plotTrajectories, ...
              Pedest_prediction, Pedest_dynamics)

% Computes the closed loop solution for the non-uniformly spaced MPC problem 
%defined by the functions
%   runningcosts:         evaluates the running costs for state and control
%                         at one sampling instant.
%                         The function returns the running costs for one
%                         sampling instant.
%          Usage: [cost] = runningcosts(x, u, x_ref,param)
%                 with  state x, reference state x_ref, control u and
%                 system parameters param
%   constraints:          computes the value of the restrictions for a
%                         sampling instance provided the data x, x_ref and u
%                         given by the optimization method. The function
%                         also considers the possible predictions of the
%                         pedestrian crossing Cross, the system parameters
%                         param, and the decision to be made by the solver
%                         to stop or accelerate stop
%                         The function returns the value of the
%                         restrictions for a sampling instance separated
%                         for inequality restrictions c and equality
%                         restrictions ceq.
%          Usage: [c,ceq] = constraints(x, x_ref, k, param, Crosses, stop)
%                 with  state x, reference state x_ref, control u,
%                 system parameters param, Predicted crossing steps
%                 Crosses, and solver decision stop
%   linearconstraints:    sets the linear constraints of the discretized
%                         optimal control problem. This is particularly
%                         useful to set control and state bounds.
%                         The function returns the required matrices for
%                         the linear inequality and equality constraints A
%                         and Aeq, the corresponding right hand sides b and
%                         beq as well as the lower and upper bound of the
%                         control.
%          Usage: [A, b, Aeq, beq, lb, ub] = linearconstraints()
%
%   system:               evaluates the linearized dynamics equation describing the
%                         process given time scale T, state vector x and control
%                         u, the linearisation state x_curr and the
%                         dynamics matrices linearisation
%                         The function returns the state vector x at the
%                         next time instant.
%          Usage: [y] = system(x, u, T, x_curr, linearisation)
%                 with time scale T, state x, control u, linearisation
%                 state x_curr and structure variable that contains the
%                 dynamics linearized matrices linearisation
% for a given number of NMPC iteration steps (mpciterations). For
% the open loop problem, the horizon is defined by the number of
% time instances N and the sampling time T. Note that the dynamic
% can also be the solution of a differential equation. Moreover, the
% initial time tmeasure, the state measurement xmeasure and a guess of
% the optimal control u0 are required.
%   Non_linear_system:  evaluates the dynamics equation without linearisation 
%                         given  state vector x and control  u.
%                         The function returns the state vector x at the
%                         next time instant.
%          Usage: x_new = Non_linear_system(x, u, param)
%                 with state x and control u, param is the system
%                 parameters
%   discret_system_matrices:  evaluates the matrices for the linearized dynamics
%                         given  state vector x and control  u and time scale T.
%                         The function returns the state vector x at the
%                         next time instant.
%          Usage: [A_d, B_d] = discret_system_matrices(x, u, T, param)
%                 with time scale T, state x and control u, param is the system
%                 parameters
%   printHeader:         Print clarifying header for selective output of closed
%                        loop data, cf. printClosedloopData
%   printClosedloopData: Print selective output of closed loop data
%   plotTrajectories:    Graphical output of the trajectories (Visualisation)
%   Pedest_prediction:  evaluates the Predicted path of the pedestrian
%                         given its initial state vector xp0 and Prediction Horizon N, 
%                         the time scale T and the Noise w.
%                         The function returns the state vector x at the
%                         next time instant.
%          Usage: Xp = Pedest_prediction(T,xp0,N,w))
%                 with time scale T, intial state vector xp0 
%                  Prediction Horizon N and the Noise w.
%   Pedest_dynamics:     evaluates the next state of the dynamics equation
%                         given its current state vector xp 
%                         the time scale T and the Noise w.
%                         The function returns the state vector y at the
%                         next time instant when noise is considered as well 
%                         as the state vector y_nmpc when no noise is
%                         considered in the pedestrian dynamics
%          Usage: [y, y_nmpc] = Pedest_dynamics(xp,T,w)
%                 with time scale T, current state vector xp
%                  and the Noise w.
%
% Arguments:
%   mpciterations:  Number of MPC iterations to be performed
%   N:              Length of optimization horizon
%   T:              Sampling interval
%   xmeasure:       State measurement of initial value of EV
%   xp_measure:     State measurement of initial value of Pedestrian
%   u0:             Initial guess of open loop control
%   x_ref:          Reference trajectory for the EV
%   param:          A structure value that stores the parameters of the
%                   Optimal Contro Problem
%   iprint:         = 0  Print closed loop data(default)
%                   = 1  Print closed loop data and errors of the
%                        optimization method
%                   = 2  Print closed loop data and errors and warnings of
%                        the optimization method
%                   >= 5 Print closed loop data and errors and warnings of
%                        the optimization method as well as graphical
%                        output of closed loop state trajectories
%                   >=10 Print closed loop data and errors and warnings of
%                        the optimization method with error and warning
%                        description
% Internal Functions:
%   applyControl:                 applies the first control element of u to
%                                 the simulated process for one sampling
%                                 interval T.
%                                 The function returns the next closed loop 
%                                 state vector xapplied
%   shiftHorizon:                 applies the shift method to the open loop
%                                 control in order to ease the restart.
%                                 The function returns a new initial guess
%                                 u0 of the control.
%   solveOptimalControlProblem:   solves the optimal control problem of the
%                                 horizon N with different time scales for the
%                                 given initial value x0, the
%                                 initial guess u0 and the predicted crossing 
%                                 steps Cross using the specified algorithm.
%                                 The function returns the computed optimal
%                                 control u, the corresponding value of the
%                                 cost function V as well as possible exit
%                                 flags and additional output of the
%                                 optimization method.
%   costfunction:                 evaluates the cost function of the
%                                 optimal control problem over the horizon
%                                 N with different time scales T and the reference path x_ref for the current
%                                 data of the optimization method x0
%                                 and u.
%                                 The function return the computed cost
%                                 function value.
%   nonlinearconstraints:         computes the value of the restrictions
%                                 for all sampling instances provided the
%                                 data x0 and u given by the
%                                 optimization method.
%                                 The function returns the value of the
%                                 restrictions for all sampling instances
%                                 separated for inequality restrictions c
%                                 and equality restrictions ceq.
%   computeOpenloopSolution:      computes the open loop solution over the
%                                 horizon N with single time scale T for the
%                                 initial values x0 as well as the
%                                 control u.
%                                 The function returns the complete open
%                                 loop solution over the requested horizon.
%   dynamics:                     evaluates the dynamic of the system for
%                                 the whole prediction Horizon N and the
%                                 various time scales
%                                 given initial valueand x0 using the control u.
%                                 The function returns all state vector x
%                                 for the whole prediction Horizon.
%   printSolution:                prints out information on the current MPC
%                                 step, in particular state and control
%                                 information as well as required computing
%                                 times and exitflags/outputs of the used
%                                 optimization method. The flow of
%                                 information can be controlled by the
%                                 variable iprint and the functions
%                                 printHeader, printClosedloopData and
%                                 plotTrajectories.
%   Pedestrian_path:              evaluates the dynamics of the Pedestrian for
%                                 the whole prediction Horizon N and the
%                                 various time scales given initial value
%                                 and xp0.
%                                 The function returns all state vector Xp
%                                 for the whole prediction Horizon.


    tol_opt = 1e-6;
    % options for fmincon
    options = optimset('Display','off',...
                'TolFun', tol_opt,...
                'MaxIter', 10000,...
                'Algorithm', 'active-set',...
                'FinDiffType', 'forward',...
                'RelLineSrchBnd', [],...
                'RelLineSrchBndDuration', 1,...
                'TolCon', 1e-6, ...
                'TolConSQP', 1e-6); 
            

    % Initilize closed loop data
    x = [];  % Position of the EV
    u = [];
    xp = []; % position of pedestrian
    cost = 0;   % running Cost
    Cost_total = cost;
    stop = [];
    
    % Applied input of the previous time step
    u_last = [0;0];

    % Get non-uniformly spaced MPC parameters
    delta_t = param.dt;
    steps = param.steps;
    w = param.noise;
    %%%%
    
    % Start of the NMPC iteration
    mpciter = 0;
    while(mpciter < mpciterations)
        % Step (1) of the NMPC algorithm:
        %   Obtain new initial value
        x0 = xmeasure;
        
        % Compute all matrices of the dynamics for different timescales
        A_d = cell(steps);       
        B_d = cell(steps);
        
        for i=1:steps
            % update timescale
            Ti = delta_t(i); 
            % update Dynamic Matrices Ai and Bi For specific time scale
            [A_i, B_i] = discret_system_matrices(x0, [0;0], Ti, param);
            % Add the new matrices to a new cell
            A_d{i} = A_i;
            B_d{i} = B_i;
        end
         linear.A = A_d;
         linear.B = B_d;
        %%%%%%
          f = Non_linear_system(x0, [0;0], param);  
          linear.f = f;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Predict Pedesterian path
        Lane = param.Lane;
        xp0 = xp_measure;
        Xp = Pedestrian_path(Pedest_prediction, N,param,xp0);    
          
        % Check which steps the pedestrian is crossing
        Crosses = zeros(N+1,1);
        for c =1:(N+1)
            if (Xp(2,c)>= Lane(1))&& (Xp(2,c)<= Lane(2)) % if Pedestrian on Crosswalk
                Crosses(c) = 1;
            end
        end
          
         %%%%%%%%%%%
    
        % Update Pedestrian position   
        xp = [ xp, xp_measure ];
        
        % Pedestrian next State
        [xp_measure,~] = Pedest_dynamics(xp0, T, w(:,mpciter+1));
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
 
        % Step (2) of the NMPC algorithm:
        %   Solve the optimal control problem
        t_Start = tic;
        [u_new, V_current, exitflag, output, stop_solution] = solveOptimalControlProblem ...
            (runningcosts, constraints, ...
             linearconstraints, system, ...
            N, x0, u0, options, x_ref(:,mpciter+1:mpciter+N+1), ...
        param , linear, u_last,Crosses,Non_linear_system);

        t_Elapsed = toc( t_Start );
        
        % compute the predicted trajectory of the EV
        x_pred = dynamics(system,N,x0,u_new, param, linear,Non_linear_system);
        
        %   Print solution
        if ( iprint >= 1 )
            printSolution(printHeader, printClosedloopData, ...
                          plotTrajectories, mpciter, x0, u_new, ...
                          iprint, ...
                          exitflag, output, t_Elapsed, param, x_ref,x_pred, xp0, Xp, Crosses);

        end


        %   Store closed loop data
        x = [ x xmeasure ]; 
        u = [ u u_new(:,1)];
        stop = [stop stop_solution];
        
        % Step (3) of the NMPC algorithm:
        %   Apply control to process
        xmeasure = applyControl(system, T, x0, u_new, linear);
        
        
       %%% Plot Velocity profile and Cost value
       cost = cost + runningcosts(xmeasure, [u_last,u_new(:,1)], x_ref(:,1), param);
       Cost_total = [Cost_total cost];
  
%        if (param.steps == 1) % Standard MPC
%            figure(3);
%             plot(0:mpciter,x(4,:),'r');
%             axis([0 mpciterations 0 13]);
%            title('Standard MPC');
%            xlabel('time instant k');
%            ylabel('velocity v (m/s)');
%           
%            figure(5)
%            plot(0:mpciter,Cost_total(1:mpciter+1),'r');
%            title('Standard MPC');
%            xlabel('time instant k');
%            ylabel('Cost function value');
%        else
%           % Non-uniformly spaced MPC 
%           figure(2);
%           plot(0:mpciter,x(4,:),'r');
%           axis([0 mpciterations 0 13]);
%           title('Non-Uniformly Spaced MPC');
%           xlabel('time instant k');
%           ylabel('velocity v (m/s)'); 
%           
%           figure(6)
%           plot(0:mpciter,Cost_total(1:mpciter+1),'r');
%           title('Non-Uniformly Spaced MPC');
%           xlabel('time instant k');
%           ylabel('Cost function value');         
%        end
%        
%        % Plot Controller Decision (Stop or Accelerate)
%        figure(4)
%        plot(0:mpciter,stop(:),'--');
%        axis([0 mpciterations -0.5 1.5]);
%        title('Does the EV choose to stop?');
%        xlabel('time instant k');
%        ylabel('Stop?');  
       
       
        %   Prepare input restart
        u0 = shiftHorizon(u_new);
        u_last = u_new(:,1);
   
        mpciter = mpciter+1  
    end
        
end

function xapplied = applyControl(system, T, x0, u, linear)
    linearisation.A = linear.A{1};
    linearisation.B = linear.B{1};
    linearisation.f = linear.f;
    xapplied = system(x0, u(:,1), T, x0, linearisation);
end


function [u, V, exitflag, output,stop_solution] = solveOptimalControlProblem ...
    (runningcosts, constraints, ...
    linearconstraints, system, N, x0, u0, ...
    options, x_ref, param, linear, u_last,Crosses,Non_linear_system)
    % Set control and linear bounds
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    for k=1:N
        [Anew, bnew, Aeqnew, beqnew, lbnew, ubnew] = ...
               linearconstraints();
        A = blkdiag(A,Anew);
        b = [b, bnew];
        Aeq = blkdiag(Aeq,Aeqnew);
        beq = [beq, beqnew];
        lb = [lb, lbnew];
        ub = [ub, ubnew];
    end

    % Solve optimization problem
    stop_solution = 0; % Variable used to store the solution the optimizer opted for
    
    stop = 0; % stop =1 car stops, stop=0 car proceeds
    if any(Crosses)
        % Pedestrian is predicted to cross
            
            % Cost for the car to pass before pedestrian crosses
            [u_acc, V_acc, exitflag_acc, output_acc] = fmincon(@(u) costfunction(runningcosts, ...
            system, N, x0, ...
            u,x_ref,param, linear, u_last,Non_linear_system), u0, A, b, Aeq, beq, lb, ...
            ub, @(u) nonlinearconstraints(constraints, ...
            system, N, x0, u, param,x_ref,linear,Crosses,stop,Non_linear_system), options);
        
            
             % Cost for the car to wait for pedestrian to cross
             stop = 1;
            [u_stop, V_stop, exitflag_stop, output_stop] = fmincon(@(u) costfunction(runningcosts, ...
            system, N, x0, ...
            u,x_ref,param, linear, u_last,Non_linear_system), u0, A, b, Aeq, beq, lb, ...
            ub, @(u) nonlinearconstraints(constraints, ...
            system, N, x0, u, param,x_ref,linear,Crosses,stop,Non_linear_system), options);

        % Decision
            if (V_stop<V_acc) % stopping the car minimizes the cost
                stop_solution = 1;
                [u, V, exitflag, output] = deal(u_stop, V_stop, exitflag_stop, output_stop);
                
            else   % accelerating the car minimizes the cost
                stop_solution = 0;
                [u, V, exitflag, output] = deal(u_acc, V_acc, exitflag_acc, output_acc);                
            end

            % If one solver returns "non feasible" solution, opt for other
            % decision
            if (exitflag_acc == -2)&&(exitflag_stop ~= -2)
                [u, V, exitflag, output] = deal(u_stop, V_stop, exitflag_stop, output_stop);
                stop_solution = 1;
            end
            if (exitflag_acc ~= -2)&&(exitflag_stop == -2)
                [u, V, exitflag, output] = deal(u_acc, V_acc, exitflag_acc, output_acc);
                stop_solution = 0;
            end
            
             
    else
        % No Pedestrianc involved (Absence of Pedestrian Constraints)
        [u, V, exitflag, output] = fmincon(@(u) costfunction(runningcosts, ...
        system, N, x0, ...
        u,x_ref,param, linear, u_last,Non_linear_system), u0, A, b, Aeq, beq, lb, ...
        ub, @(u) nonlinearconstraints(constraints, ...
        system, N, x0, u, param,x_ref,linear,Crosses,stop,Non_linear_system), options);
        
    end

end

function x = computeOpenloopSolution(system, N, T, x0, u, linearisation)
    x(:,1) = x0;
    for k=1:N
        x(:,k+1) = system(x(:,k), u(:,k), T, x0, linearisation);                       
    end
end

function u0 = shiftHorizon(u)
    u0 = [u(:,2:size(u,2)) u(:,size(u,2))];
end

function cost = costfunction(runningcosts, system, ...
                    N, x0, u, x_ref, param, ...
                    linear, u_last,Non_linear_system)
    cost = 0;
    n = length(x0);
    x = zeros(n, N+1);
    x = dynamics(system,N,x0,u, param, linear,Non_linear_system); 
    cost = cost+ runningcosts(x(:,2), [u_last,u(:,1)], x_ref(:,2), param);
    
    for k=2:N
        cost = cost+runningcosts(x(:,k+1), u(:,k-1:k), x_ref(:,k+1), param);
    end
end

function [c,ceq] = nonlinearconstraints(constraints, ...
    system, ...
    N, x0, u, param,x_ref,linear,Crosses,stop, Non_linear_system)
    x = zeros(N+1, length(x0));

    x = dynamics(system,N,x0,u, param, linear, Non_linear_system);

    c = [];
    ceq = [];
    for k=2:N
        [cnew, ceqnew] = constraints(x(:,k),x_ref(:,k), k, param,Crosses,stop);
        c = [c cnew];
        ceq = [ceq ceqnew];
    end
    % terminal constraints
    [cnew, ceqnew] = constraints(x(:,N+1),x_ref(:,N+1), N+1, param,Crosses,stop);
    c = [c cnew];
    ceq = [ceq ceqnew];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% OUTPUT %%%%%%%%%%%%%
function printSolution(printHeader, printClosedloopData, ...
             plotTrajectories, mpciter, x0, u, iprint, exitflag, output, t_Elapsed, ...
             param, x_ref, x, xp0, Xp, Crosses)
    % Print results
    if (mpciter == 0)
        printHeader();
    end
    printClosedloopData(mpciter, u, x0, t_Elapsed);
    % Print errors and warnings (fmincon)
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
    
        % Plot the trajectory
        plotTrajectories(x, x0, param, x_ref, xp0, Xp, Crosses) 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  NON UNIFORM NPMC  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = dynamics(system, N, x0, u0, param, linear,Non_linear_system)
%% Calculate the state vector X 
% Compute the state vector x by Concatenating all state vectors
% corresponding to different timescales

    %parameters
    A_d = linear.A;    % Matrices A for all timescales (bicycle model Dynamics)
    B_d = linear.B;    % Matrices B for all timescales (bicycle model Dynamics)
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
    
        % Calculate timescale corresponding state vector
        xi = computeOpenloopSolution(system, Np, T, x0, u0(:,1+Np*(i-1):Np*i), linearisation);
        % update parameters
        x0 = xi(:,Np+1); 
        xi = xi(:,2:Np+1); % remove x(:,1)=x0
        
        linearisation.f = Non_linear_system(x0, [0;0], param); 

        % Concatenate new vectors
        X = [X,xi];
     
end

end

function Xp = Pedestrian_path(Pedest_prediction, N,param,xp0)
%% Calculate the non uniform pedestrian state Xp 
% Compute the state vector x by Concatenating all state vectors
% corresponding to different timescales
    delta_t = param.dt;
    steps = param.steps;
    w = param.noise;
    Xp = xp0;
    Np = N/steps; %length of prediction timescale
    for i=1:steps
        % update timescale
        T = delta_t(i);
        
        % Calculate timescale corresponding state vector
        xpi = Pedest_prediction(T,xp0,Np,w);
        
        % update parameters
        xp0 = xpi(:,Np+1); 
        xpi = xpi(:,2:Np+1); % remove xp(:,1)=xp0
        
        % Concatenate new vectors
        Xp = [Xp,xpi];        
    end

    

end