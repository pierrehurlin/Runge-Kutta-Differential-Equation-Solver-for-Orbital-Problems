%------------------------------------------------------------------------------
% PURPOSE: This code reproduces Figure 1 and produces four ellipses with
% different eccentricity values  
% -----------------------------------------------------------------------------
% Author: Pierre Hurlin 
% DKE, Maastricht University
% Version v1.01, July 26, 2021 
% Resit - Project 1.2
% ---------------------------------------------------------------------

clear ; clc ; close all 

%==================
%=== Parameters ===
%==================
a=1;                        % Semi-major axis
E=[0 0.25 0.5 0.85] ;       % Eccentricity values

figure 

for i=1:4
    
    e=E(i);                     % Eccentricity
    b = a*sqrt(1-e^2) ;         % Semi-minor axis 
    c=sqrt(a^2-b^2);            % Distance from the center to the focus S
    q=a*(1-e);                  % Distance from the focus to the perikorn

    subplot(2,2,i)
    C = [-c 0] ;                 % Center 
    th = linspace(0,2*pi);  
    xe = C(1)+a*cos(th) ; 
    ye = C(2)+b*sin(th) ; 
    plot(xe,ye,'b','Linewidth',1.5)
    title(sprintf('Eccentricity = %1.2f ',e))
    grid('on')
    hold('on')
    ylim([-1.5*b 1.5*b])
    axis equal
    plot(0,0,'o','Color','r','Linewidth',1.5)
    %legend('Ellipse','Focus S','Location','best')
    %legend('boxoff')
end
