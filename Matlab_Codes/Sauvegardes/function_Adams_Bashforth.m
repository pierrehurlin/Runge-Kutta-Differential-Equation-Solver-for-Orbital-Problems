%------------------------------------------------------------------------------
% PURPOSE: solve a differential equation using the 3rd-order Adams-Bashforth method 
% with a fourth-order Runga-Kutta boostrap.
% -----------------------------------------------------------------------------
% function [omega,omega_RK,k]=function_Adams_Bashforth(h,y0,t,f)
%
% INPUTS:
%   h: step-size 
%   y0: initial condition on y(t) 
%   t: vector of time 
%   f: function y_point(t)
%
% OUTPUT: 
%   omega is a Tx1 vector of values y(t) obtained with the 3rd-order Adams-Bashforth method
%   omega_RK is a 3x1 vector of values y(t) obtained with a fourth-order Runge-Kutta method
%   k: 3x4 matrix of components of the fourth-order Runge-Kutta method
%
% ---------------------------------------------------------------------- 
% P. Hurlin. April 22, 2021.
% DKE, Maastricht University
% ---------------------------------------------------------------------

function [omega,omega_RK,k]=function_Adams_Bashforth(h,y0,t,f)

T=length(t);                % Number of periods

%========================================
%=== Fourth-order Runge-Kutta methods ===
%========================================
% Comment: for efficiency, we only compute the Runge-Kutta values necessary 
% to bootsrap the three-stage Adams-Bashforth method. 
% Otherwise please change the T_RK value from 3 to T
%
T_RK=3;                                    % Number of initial conditions necessary for the Bootsrap (or T_RK=T for all values)
k1=NaN(T_RK,1);k2=k1;k3=k1;k4=k1;          % Initialization of the k1, k2, k3 and k4 components
omega_RK=NaN(T_RK,1);omega_RK(1)=y0;       % Initialization of the omega vector 

for i=1:T_RK-1
    k1(i)=h*f(t(i),omega_RK(i));
    k2(i)=h*f(t(i)+0.5*h,omega_RK(i)+0.5*k1(i));
    k3(i)=h*f(t(i)+0.5*h,omega_RK(i)+0.5*k2(i));
    k4(i)=h*f(t(i)+h,omega_RK(i)+k3(i));
    omega_RK(i+1)=omega_RK(i)+(1/6)*(k1(i)+2*k2(i)+2*k3(i)+k4(i));
end
k=[k1 k2 k3 k4];

%========================================
%=== 3rd-order Adams-Bashforth method ===
%========================================
omega=NaN(T,1);                    % Initilization of the omaga vector
omega(1:3)=omega_RK(1:3);          % Boostrapping with the Runge-Kutta values
for i=3:T-1
    omega(i+1)=omega(i)+(1/12)*h*(23*f(t(i),omega(i))-16*f(t(i-1),omega(i-1))+5*f(t(i-2),omega(i-2)));
end
