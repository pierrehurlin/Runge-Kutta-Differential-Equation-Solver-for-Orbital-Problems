%------------------------------------------------------------------------------
% PURPOSE: Solve a differential equation using the 3rd-order Adams-Bashforth method 
% with a fourth-order Runga-Kutta boostrap and a step-size of h=0.1 and h=0.05
% Use the secant method to find the time at which y(t) = 0 to an accuracy of 1e-03.
% -----------------------------------------------------------------------------
% Author: Pierre Hurlin 
% DKE, Maastricht University
% Version v1.01, April 22, 2021 
% Numerical Mathematics DKE - Assignement  
% ---------------------------------------------------------------------
% Note: here, I use a fourth-order Runge-Kutta method as it is recommended 
% to use a method of the same order (or higher) to boostrap the multistep method, 
% i.e., the 3rd-order Adams-Bashforth.
% ---------------------------------------------------------------------

clear ; clc ; close all ;
help Assignment_Pierre_Hurlin_DKE

%================
%== Parameters ==
%================
t0=0;tf=6;                  % Initial and final time      
y0=2;                       % Initial value!!
f=@(t,y) sin(t)+y-y.^3;     % y_point(t) 
True_value=-0.6599693029;   % True value

fprintf('   Initial value: y(0) = %1.0f \n',y0)
fprintf('   Time span: t0 = %1.0f  tf = %1.0f \n',t0,tf )

%===========================
%== Loop on the step-size ==
%===========================

for h=[0.1 0.05]                            % Loop on the step-size

    t=(t0:h:tf)';                           % Time
    T=length(t);                            % Number of periods
    ii=(0:1:T-1)';                          % Index

    % Comment: the function_Adams_Bashforth computes the 3rd-order Adams-Bashforth 
    % method with a fourth-order Runge-Kutta boostrap

    [omega,omega_RK,k]=function_Adams_Bashforth(h,y0,t,f);      

    disp(' '),fprintf('   ==== Step-size: h = %1.2f ==== \n',h)

    disp(' '),disp('    ---- Boostrap: Fourth-order Runge-Kutta method ----'), disp(' ')
    varNames = {'i','t','omega','k1','k2','k3','k4','k'};
    Runge_Kutta=table((0:1:size(k,1)-1)',t(1:size(k,1)),omega_RK,k(:,1),k(:,2),k(:,3),k(:,4),sum(k,2),'VariableNames',varNames);
    disp(Runge_Kutta)

    disp(' '),disp('    ---- Three-stage Adams-Bashforth method ----'), disp(' ')
    varNames = {'i','t','omega'};
    Adams_Bashforth=table(ii((t<=0.5)|(t==tf)),t((t<=0.5)|(t==tf)),omega((t<=0.5)|(t==tf)),'VariableNames',varNames);
    disp(Adams_Bashforth)

    absolute_error=abs(True_value-omega(end));
    disp('    ---- Absolute error ---- ')
    fprintf('   True value: y(6) = %1.10f \n',True_value)
    fprintf('   Approximate value (Three-stage Adams-Bashforth): y(6) = %1.10f \n',omega(end))
    fprintf('   Absolute error = %1.10f \n',absolute_error)
    fprintf('   Relative error = %1.4f %% \n',abs(absolute_error/True_value)*100)

end

%======================================
%=== Find the time at which y(t)=0 ==== 
%======================================
omega_inf=1;                                            % Initialisation

while abs(omega_inf)>1e-03                              % Loop on the accuracy
    
    t_inf=t(find(omega>0,1,'last'));                    % Last time t with a positive y(t)
    t_sup=t(find(omega<0,1,'first'));                   % First time t with a negative y(t)
    omega_inf=omega(find(omega>0,1,'last'));            % y(t_inf)
    omega_sup=omega(find(omega<0,1,'first'));           % y(t_sup)
    
    t_cut=t_sup-(t_sup-t_inf)/(omega_sup-omega_inf)*omega_sup;  % Secant method
    h_k=t_cut-t_inf;                                            % Adjusted step-size
    t=(t0:h_k:tf)';                                             % New grid
   
    [omega,omega_RK,k]=function_Adams_Bashforth(h_k,y0,t,f);    %  3rd-order Adams-Bashforth method

end

disp(' '), disp('   ===== Find the time at which y(t)=0 ===')
fprintf('   t = %1.5f \n',t_cut)
fprintf('   h_k = %1.5f \n',h_k)
fprintf('   y(t_inf) = %1.5f \n',omega_inf)
fprintf('   y(t_sup) = %1.5f \n',omega_sup)
