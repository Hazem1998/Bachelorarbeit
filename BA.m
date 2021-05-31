function Main
%% Clear and addpath
    clc;
    close all;
    clear all;
    addpath functions;
%% General
    mpciterations = 40;     % Simulation length
    N             = 12;     % prediction Horizon  
    T             = 0.1;    % Sampling interval
    iprint = 11;            % Print and visualize Output
    
%% Non uniform parameters
    steps = 3; % number of timescales (non uniform steps)
    
    delta_t = size(steps,1);
    delta_t(1) = T;
    for i=2:steps
        delta_t(i) = delta_t(i-1) + T; % Time scales (as per each non-uniformly spaced MPC step)
    end
    param.dt = delta_t;
    param.steps = steps;
    
    
%% Model parameters
    l_r = 1.5;    % Distance from vehicle center of gravity to the rear Axis in meters  
    l_f = 1.5;    % Distance from vehicle center of gravity to the front Axis in meters
    param.distance = [l_r, l_f];
    % Weighting matrices (Cost function)
    param.Q = diag([0 15 15 0.2]);
    param.R = diag([0.33 15]);
    param.S = diag([5 5]);

    % Lane coordinates 6.5m width
    param.Lane = [-1.5,5]; 

    % staying in Lane constraints
    param.dev = 1; % distance the car is allowed to diviate from the reference trajectory

    % Safety parameters in case of pedestrian crossing
    param.crossing = [18,20]; % pedestrian crossing coordinates
    param.safety = 1; % distance of pedestrian from lane to consider stopping
    param.v_max = 13; % the max velocity of the car 15ms = 50km/h

%% Initializations
    xmeasure      = [0.0; 0.0; 0.0; 8.5];  % Initial state
    u0            = zeros(2,N);  % Input initial guess
    xp_measure = [19;-5;0;1.5];    % initial position of pedestrian
    
%% reference trajectory
    x_ref = [0:mpciterations+N;zeros(2,mpciterations+N+1);8.5*ones(1,mpciterations+N+1)];

%% Optimization
    % Define same noise for both simulations
    param.noise = randn(2,mpciterations);

    % Non uniform MPC
    nmpc_me(@runningcosts, @constraints, ...
          @linearconstraints, @system, ...
         mpciterations, N, T, xmeasure, u0, x_ref,xp_measure,iprint,@Non_linear_system, ...
            param, @discret_system_matrices, @printHeader, @printClosedloopData, @plotTrajectories, ...
        @Pedest_prediction, @Pedest_dynamics);
    
    % Standard MPC
    param.steps = 1;
    
    nmpc_me(@runningcosts, @constraints, ...
          @linearconstraints, @system, ...
         mpciterations, N, T, xmeasure, u0, x_ref,xp_measure,iprint,@Non_linear_system, ...
            param, @discret_system_matrices, @printHeader, @printClosedloopData, @plotTrajectories, ...
        @Pedest_prediction, @Pedest_dynamics);

    rmpath('./functions');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the NMPC functions
% x=[s;d;phi;v] size(x) = 4 x 1
% u=[a;delta]   size(u) = 2 x N with N= Horizon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(x, u, x_ref,param)
%% Running cost of the system
   u_curr = u(:,2);
   u_prev = u(:,1);
   Q = param.Q;
   S = param.S;
   R = param.R;
   cost = (x-x_ref).'*Q*(x-x_ref) + u_curr.'*R*(u_curr) + (u_curr-u_prev).'*S*(u_curr-u_prev);
end


function [c,ceq] = constraints(x, x_ref, k, param, Crosses, stop)
%% Non linear constraints (State Constraints)
    % Staying in lane conditions
    w_lane = param.dev; 
    L = abs( x(2)-x_ref(2) )-w_lane;
    
    % velocity conditions
    vmax = x(4)- param.v_max;
    v_positive = -x(4);
    
    % physical constraints
    xp_lim = param.crossing;
    l_r = param.distance(1);
    l_f = param.distance(2);
    
    % Constraints vector: 
    % Pedestrian Constraint == 0 unless a pedestrian is detected
    c   = [L;vmax;v_positive;-1];
    
    % The pedestrian is predicted to cross
    if any(Crosses) 

        %if the car is stopping before the crossing
        if (stop ==1)  
            % Find the furthest predicted step the pedestrian is crossing
            i_step = find(Crosses,1,'last');
            % apply the Stop constraints for all steps before that
            if (k<=i_step) 
                c(end) = x(1)- (xp_lim(1)-l_f);
            end           
        end
        
        %if the car is accelerating
        if (stop ==0)
            % Apply the Acceleration Constraint ONLY for the steps the
            % pedestrian is crossing
            if (Crosses(k) == 1)
                c(end) = -x(1)+ (xp_lim(2)+l_r);
            end     
        end
    end
    
    ceq = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints()
%% Input Constraints
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = [-9; -0.2];   % Min input
    ub  = [5;0.2];      % Max input
end



function x_new = Non_linear_system(x, u, param)
 %% Dynamics of the non linearized system: update current state but not for prediction
 [~,~,phi,v,a,delta] = curr_state_input(x, u);
 
 l_r = param.distance(1);
 l_f = param.distance(2);
 
 alpha = atan( (l_r* tan(delta) ) / (l_f + l_r) );
 s_new = v*cos(phi + alpha);
 d_new = v*sin(phi + alpha);
 phi_new = (v*sin(alpha)) / l_r;
 v_new = a;

 x_new = [s_new;d_new;phi_new;v_new];


end


 function y = system(x, u, T, x_curr, linearisation)
 %% Dynamics of the system 
    A_d = linearisation.A;
    B_d = linearisation.B;
    f = linearisation.f;

    y = x_curr + T*f + A_d*(x - x_curr) + B_d*u;
 end
 
 function [A_d, B_d] = discret_system_matrices(x, u, T, param)
 %% Dynamic system Matrices
 [~,~,phi,v,~,delta] = curr_state_input(x, u);

 l_r = param.distance(1);
 l_f = param.distance(2);
 
 
 z1 = phi + atan(l_r*tan(delta) / (l_r + l_f));
 z2 = T*T*v*tan(delta); 
 z3 = ( (l_r*tan(delta)) / (l_r+l_f) )^2  + 1;
 z4 = (l_r + l_f)*sqrt(z3);
 z5 = (l_r + l_f)^3 * z3^(3/2);
 z6 = (l_r + l_f)*z3;
 z7 = v*(1 + tan(delta)^2 )*(1/z4 - ( (l_r*tan(delta))^2 / z5) ) ;
 z8 = T*l_r*v*(1 + tan(delta)^2 );
 
 
 A13 = -T*v*sin(z1);
 A14 = T*cos(z1) - ( z2*sin(z1)/(2*z4) );
 A23 = T*v*cos(z1);
 A24 = T*sin(z1) + ( z2*cos(z1)/(2*z4) );
 A34 = T*tan(delta)/z4;
 

 
 B = zeros(4,2);
 B(1,1) = T*T*cos(z1)/2;
 B(1,2) = -T*T*v*z7*sin(z1)/2 - (z8*sin(z1)/z6);
 B(2,1) = T*T*sin(z1)/2;
 B(2,2) = T*T*v*z7*cos(z1)/2 + (z8*cos(z1)/z6);
 B(3,1) = T*T*tan(delta)/(2*z4);
 B(3,2) = T*z7;
 B(4,1) = T;
 B(4,2) = 0;
 
 A_d = [1 0 A13 A14; 0 1 A23 A24; 0 0 1 A34; 0 0 0 1];
 B_d = B;
 
 end
    
 function [s,d,phi,v,a,delta] =  curr_state_input (x, u)
 %% Split the state and input vectors into components
 s = x(1);
 d = x(2);
 phi = x(3);
 v = x(4);
 a = u(1);
 delta = u(2);
 end

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% OUTPUT %%
 function printHeader()
    fprintf('   k  |      a(k)        delta(k)        s        d        phi        v     Time\n');
    fprintf('--------------------------------------------------------------------------------\n');
end

function printClosedloopData(mpciter, u, x, t_Elapsed)
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+6.3f', ...
             mpciter, u(1,1),u(2,1), x(1), x(2),x(3),x(4), t_Elapsed);
end

function plotTrajectories( x, x0, param, x_ref, xp0, Xp, Crosses)           
%% Plot current frame
    l_r = param.distance(1);
    l_f = param.distance(2);
    L = l_r + l_f;  % Vehicule length
    B = L/2;        % vehicule width
    
    C = [0. 0. ; L 0. ; L B ; 0. B; 0 0.] ; % vehicule coordinates
    %%%%%%%%%%%%%%%%%%%%%%%
    xp_lim = param.crossing; 
    y_lim = param.Lane;
    y_middle = (y_lim(1)+y_lim(2))/2; % middle line between lanes
    
    % Color the pedestrian path
    x_ped = [xp_lim(1), xp_lim(1), xp_lim(2), xp_lim(2), xp_lim(1)];
    y_ped = [y_lim(1), y_lim(2), y_lim(2), y_lim(1), y_lim(1)];
    %%%%%%%%%%%
    figure(1);
        title('Velocity v = 10 (m/s)');
        xlabel('position x (m)');
        ylabel('position y (m)');
        grid on;
        hold on;
        
        %Plot pedestrian path
        patch(x_ped, y_ped, 'black', 'FaceColor', 'black', 'FaceAlpha', 0.1);
       
        % Plot Lane
        plot(linspace(-2,30,100),y_lim(1)*ones(100),'black')
        plot(linspace(-2,30,100),y_middle*ones(100), '--')
        plot(linspace(-2,30,100),y_lim(2)*ones(100),'black')
        
        
        % Plot Pedesterian predictions
        for i=1:size(Crosses)
            if (Crosses(i) == 1)
                plot(Xp(1,i),Xp(2,i),'xr', 'MarkerFaceColor','g')
                plot(x(1,i),x(2,i),'or', 'MarkerFaceColor','r');
            else
                plot(Xp(1,i),Xp(2,i),'xb', 'MarkerFaceColor','g')
                plot(x(1,i),x(2,i),'ob');
            end
        end
        
        % Plot Pedesterian current position
        plot(xp0(1),xp0(2),'o','MarkerFaceColor','black')
        
        % Plot reference Path
        plot(linspace(x_ref(1,1),x_ref(1,length(x_ref)),100),linspace(x_ref(2,1),x_ref(2,length(x_ref)), 100),'r')
              
       % Plot the actual trajectory of the car
        plot(C(:,1)+x0(1)-L/2,C(:,2)+x0(2)-B/2,'b')
        
        axis([-2 30 -5 7]);
        pause(1)
        close   % used to update each frame      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Pedestrian%%%%%%%
% xp = [x,y,vx,vy] 4 x 1 
function [y, y_nmpc] = Pedest_dynamics(xp,T,w) 
%% Pedestrian dynamics

    F = [1 0 T 0; 0 1 0 T; 0 0 1 0; 0 0 0 1];
    t = (T^2)/2;
    G = [0 0;0 t; 0 0; 0 T];
    
    % Prediction of the pedestrian next movement
    y_nmpc = F*xp;

    % Actual movement with Gaussian Noise
    %w = randn(2,1);
    y = y_nmpc + G*w;
end

function Xp = Pedest_prediction(T,xp0,N,w)
%% Prediction of the pedestrian Path
    Xp(:,1) = xp0;
    for k=1:N
       [~,yp]  = Pedest_dynamics(Xp(:,k),T,w(:,1));  
        Xp(:,k+1)= yp;   
    end
end









