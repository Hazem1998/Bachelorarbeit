function [tangent,normal,k,x,y] = Bez_current(s, P, n)
% function returns the parameters of Bezier curve at a given s
% P=[P1; P2;...; Pn] Control points with P1 = [x1 y1]
% n: number of points , s: curvi abscissa
t_s = inv_abscissa(s,t,P,n);
Bez = Bezier(t_s, P, n) ;
[x,y] = Bez(t_s,:);
% tangent
Bez_deriv = deriv_Bezier(t_s,P,n);
tangent = Bez_deriv(t_s,:);
tangent = tangent/norm(tangent,2) % normalaize the tangent
% calculate the normal
normal = (-tangent(2),tangent(1));
% Curvature
%k = how to get deriv of x and y 
% dt = 0.001
% x_last = Bez(t_s - dt,:)
% x_der = (x - x_last)/dt

end




%%%%%%%%%%%%%%%%
%%%%% Bezier%%%%
%%%%%%%%%%%%%%%%

function Curve = Bezier(t, P, n) 
% Function returns the Bezier Curve given control points
% t: time
for    i=0:1:(n-1)
binomial(i+1)=factorial(n-1)/(factorial(i)*factorial(n-1-i));  % for calculating (x!/(y!(x-y)!)) values 
end

l=[];
B_i=[];
for tj=0:0.001:t % tj: samples from 0 to t
for j=1:n
B_i(j)=binomial(j)*((1-tj)^(n-j))*(tj^(j-1)); % Bezier curve formel B_i
end
l=cat(1,l,B_i);  %catenation: % size[t/0.002 x n]
end
Curve=l*P;  % size[t/0.002 x 2] == [t ==> x,y]

end

function Curv_deriv = deriv_Bezier(t,P,n)
% the derivative of Bezier curve is Bezier curve n-1 points
Q = [];
for i =1:(n-1)
  Q(i) = n(P(i+1) - P(i)); 
  %Q(i) = n(P(i+1) - P(i))/ 0.001; % Change I found 2 formels online that is why
end
Curv_deriv = Bezier(t, Q, n-1);

% Curve = Bezier(t1, P, n); % Bezier Curve 
% %[t1,~]=size(Curve);
% t0 = t1 -0.001;
% p1 = Curve(t1,:);
% p0 = Curve(t0,:);
% Curv_deriv = n*(p1 - p0)/(t1 - t0);
% if t1 = 0 % Change: i think it is not possible to start from 0 since there will be no derivative,
%     %maybe calculate the limit by hand
%     Curv_deriv = 
% end

end

function s = abscissa(t,P,n)
% Calculate abscissa
s =0;
Curv_deriv = deriv_Bezier(t,P,n);
for tau = 0:0.001:t 
 s = s + norm(Curv_deriv(tau,:),2);
end

end

function t_s = inv_abscissa(s,t,P,n)
% Compute inverse of abscissa t(s) from the function s = abscissa(t,P,n)
% still not able to figured out how

end
