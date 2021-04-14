function Main
%% Clear and addpath
clc;
close all;
clear all;
addpath functions;
%% General
    mpciterations = 5;     % How this was chosen?
    N             = 20;     
    T             = 0.1;    % Sampling interval
    
%% Model parameters
l_r = 1.4;    % Distance from vehicle center of gravity to the rear in meters
l_f = 1.4;    % Distance from vehicle center of gravity to the front in meters
u_init = [0;0]; % from PAPER: starting with a non zero input results large deviation for far prediction
param.distance = [l_r, l_f];
param.input = u_init;

%% Initializations
    tmeasure      = 0.0;
    xmeasure      = [0.0, 0.0, 0.0, 0.0];  % starts from equilibrium
    u0            = ones(2,N);  % this is initial guess
%% reference trajectory
x_ref = [transpose(1:mpciterations),zeros(mpciterations,2),5*ones(mpciterations,1)];

%% Optimization

    nmpc_me(@runningcosts, @constraints, ...
         @terminalconstraints, @linearconstraints, @system, ...
         mpciterations, N, T, tmeasure, xmeasure, u0, x_ref,@Non_linear_system, ...
            param, @discret_system_matrices, @printHeader, @printClosedloopData, @plotTrajectories);

    rmpath('./functions');

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Definition of the NMPC functions
% x=[s,d,phi,v] size(x) = 1 x 4 
% u=[a;delta]   size(u) = 2 x N with N= Horizon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x, u, x_ref)
%% Running cost of the system
   u_curr = u(:,2);
   u_prev = u(:,1);
    cost = norm(x-x_ref,2)+norm(u,2) + norm(u_curr-u_prev,2); % CHANGE = %dont forget to add the Q,R and S weight matrices (are they given as Input? / initialization?) //

end


function [c,ceq] = constraints(t, x, u)
    c   = [];
    ceq = [];
end

function [c,ceq] = terminalconstraints(t, x)
    c   = [];
    ceq = [];
end

function [A, b, Aeq, beq, lb, ub] = linearconstraints(t, x, u)
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = [];
    ub  = []; 
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

x_new = [s_new,d_new,phi_new,v_new];


end


 function y = system(t, x, u, T, x_curr, param, linearisation)
 %% Dynamics of the system 
u_init = param.input;
A_d = linearisation.A;
B_d = linearisation.B;
f = linearisation.f;

y = x_curr + T*f + transpose(A_d*(x' - x_curr') + B_d*(u - u_init));
 end
 
 function [A_d, B_d] = discret_system_matrices(t, x, u, T, param)
 %% Values taken from "Appendix A", the LINEARIZED AND DISCRETIZED SYSTEM MATRICES
 [~,~,phi,v,~,delta] = curr_state_input(x, u);

 l_r = param.distance(1);
 l_f = param.distance(2);
 
 
 z1 = phi + atan(l_r*tan(delta) / (l_r + l_f));
 z2 = T*T*v*tan(delta); 
 z3 = (l_r*tan(delta))^2 / (l_r+l_f)^2;
 z4 = (l_r + l_f)*sqrt(z3);
 z5 = (l_r + l_f)^3 * z3^(3/2);
 z6 = (l_r + l_f)*z3;
 z7 = v*(1 + tan(delta)^2 )*(1/z4 + (l_r*tan(delta)^2/z5) ) ;
 z8 = T*l_r*v*(1 + tan(delta)^2 );
 
 
 A13 = -T*v*sin(z1);
 A14 = T*cos(z1) - z2*sin(z1)/2*z4;
 A23 = T*v*cos(z1);
 A24 = T*sin(z1) + z2*cos(z1)/2*z4;
 A34 = T*tan(delta)/z4;
 

 
 B = zeros(4,2);
 B(1,1) = T*T*cos(z1)/2;
 B(1,2) = -T*T*v*z7*sin(z1)/2 - z8*sin(z1)/z6;
 B(2,1) = T*T*sin(z1)/2;
 B(2,2) = T*T*v*z7*cos(z1)/2 + z8*cos(z1)/z6;
 B(3,1) = T*T*tan(delta)/(2*z4);
 B(3,2) = T*z7;
 B(4,1) = T;
 B(4,2) = 0;
 
 A_d = [1 0 A13 A14; 0 1 A23 A24; 0 0 1 A34; 0 0 0 1];
 B_d = B;
 
 end
    
 function [s,d,phi,v,a,delta] =  curr_state_input (x, u )
 %% Split the state and input vectors into components
 s = x(1);
 d = x(2);
 phi = x(3);
 v = x(4);
 a = u(1);
 delta = u(2);
 end

 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% OUTOUT %%
 function printHeader()
    fprintf('   k  |      a(k)        delta(k)        s        d        phi        v     Time\n');
    fprintf('--------------------------------------------------------------------------------\n');
end

function printClosedloopData(mpciter, u, x, t_Elapsed)
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+11.6f %+6.3f', ...
             mpciter, u(1,1),u(2,1), x(1), x(2),x(3),x(4), t_Elapsed);
end

function plotTrajectories(system, T, t0, x0, u, param, linearisation, x_ref)
                          
    x = system(t0, x0, u, T, x0, param, linearisation);
    
    figure(1);
        title('s/d closed loop trajectory');
        xlabel('s');
        ylabel('d');
        grid on;
        hold on;
        plot(x_ref(:,1),x_ref(:,2),'o', ...
        linspace(x_ref(1,1),x_ref(length(x_ref),1),100),linspace(x_ref(1,2),x_ref(length(x_ref),2), 100),'r')
        plot(x(1),x(2),'or', ...
             'MarkerFaceColor','r');
        axis([-2 6 -3 3]);
        axis square;

end
