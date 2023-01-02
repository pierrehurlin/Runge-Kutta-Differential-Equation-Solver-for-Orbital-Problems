%------------------------------------------------------------------------------
% PURPOSE: This code implements a 3rd-order Runge-Kutta differential equation solver 
% and applies it your solver to the problem of a satellite orbiting a large central body 
% -----------------------------------------------------------------------------
% Author: Pierre Hurlin 
% DKE, Maastricht University
% Version v1.03, August 3, 2021 
% Resit - Project 1.2
% ---------------------------------------------------------------------
% Note: we choose here the Heun's third order method
% ---------------------------------------------------------------------

clear ; clc ; close all 
help Code_Resit_Project_12

%==================
%=== Parameters ===
%==================
mu=1;                       % Product of the gravitational constant by the mass of Saturn 
a=1;                        % Semi-major axis
T=2*pi*sqrt((a^3)/mu);      % Period 
e=0.5;                      % Eccentricity
b = a*sqrt(1-e^2) ;         % Semi-minor axis 
c=sqrt(a^2-b^2);            % Distance from the center to the focus S
q=a*(1-e);                  % Distance from the focus to the perikorn

t0=0;                       % Initial time
tf=T;                       % Final time
h=0.1;                      % Step size
n=floor((tf-t0)/h)+1;       % Number of steps
t=(t0:h:t0+n*h)';           % Time
  
fprintf('   =================== \n')
fprintf('   ==== Parameters === \n')
fprintf('   =================== \n \n')
fprintf('   Product GxM = mu = %1.4f \n',mu)
fprintf('   Semi-major axis a = %1.4f \n',a)
fprintf('   Semi-minor axis b = %1.4f \n',b)
fprintf('   Distance focus S c = %1.4f \n',c)
fprintf('   Eccentricity e = %1.4f \n',e)
fprintf('   Period T = %1.4f \n \n',T)

fprintf('   Step size h = %1.4f \n',h)
fprintf('   Number of steps required to et a full period n = %1.0f \n',n)
fprintf('   Time span: t0 = %1.0f  tf = %1.4f \n',t0,tf )

%==========================
%=== Initial conditions ===
%==========================
x0=q;
y0=0;
d=sqrt(x0^2+y0^2); 
vx0=0;
vy0=sqrt(mu/a*(1+e)/(1-e));

z0=[x0 y0 vx0 vy0];                       % Vector of initial conditions
disp(' ')
fprintf('   =========================== \n')
fprintf('   ==== Initial conditions === \n')
fprintf('   =========================== \n \n')
fprintf('   x_{0} = %1.0f \n',z0(1))
fprintf('   y_{0} = %1.0f \n',z0(2))
fprintf('   v_{x,0} = %1.0f \n',z0(3))
fprintf('   v_{y,0} = %1.0f \n',z0(4))

%=============================================
%=== Newton's law of universal gravitation ===
%=============================================
f=@(t,z) [z(3)  z(4)  -mu/((z(1)^2+z(2)^2)^(3/2))*z(1)  -mu/((z(1)^2+z(2)^2)^(3/2))*z(2)] ;     % f(t,z) 

%=============================
%=== RK third order method ===
%=============================
% Note: here, we consider the Heun's method, but it is possible to choose
% other c2 and c3 parameters to get automatically any other RK 3rd order method

k1=NaN(n+1,4);k2=k1;k3=k1;                 % Initialization of the k1, k2, and k3 components
w_RK=NaN(n+1,4);w_RK(1,:)=z0;              % Initialization of the w vector 
c2=1/3;c3=2/3;                             % Parameters for Heun's method (change for any other 3rd RK method)

b3=(3*c2-2)/(6*c3*(c2-c3));                % Constrained parameters for 3rd RK method
b2=(3*c3-2)/(6*c2*(c3-c2));                % Constrained parameters for 3rd RK method
b1=1-b2-b3;                                % Constrained parameters for 3rd RK method
a32=1/(6*c2*b3);a31=c3-a32;a21=c2;         % Constrained parameters for 3rd RK method

disp(' ')
fprintf('   =============================== \n')
fprintf('   ==== Runge Kutta Parameters === \n')
fprintf('   =============================== \n \n')
fprintf('   c2 = %1.4f  c3 = %1.4f  \n',c2 , c3)
fprintf('   b1 = %1.4f  b2 = %1.4f  b3 = %1.4f \n',b1,b2,b3)
fprintf('   a21 = %1.4f  a31 = %1.4f  a32 = %1.4f \n',a21,a31,a32)

w=NaN(n+1,4);w(1,:)=z0;                % Initialization of the w vector (naive)

for i=1:n
    k1(i,:)=h*f(t(i),w_RK(i,:));
    k2(i,:)=h*f(t(i)+c2*h,w_RK(i,:)+a21*k1(i,:));
    k3(i,:)=h*f(t(i)+c3*h,w_RK(i,:)+a31*k1(i,:)+a32*k2(i,:));
    w_RK(i+1,:)=w_RK(i,:)+(b1*k1(i,:)+b2*k2(i,:)+b3*k3(i,:));
end

res=[t w_RK];                           % Results
x=res(:,2);                             % x coordinate
y=res(:,3);                             % y coordinate
vx=res(:,4);                            % vx velocity x dimension 
vy=res(:,5);                            % vy velocity y dimension 

disp(' ')
fprintf('   ============================ \n')
fprintf('   ==== Runge Kutta results === \n')
fprintf('   ============================ \n \n')

if n<=10                      % Display parameter: we report the dis first lines of the result matrix
    dis=n+1; 
else
    dis=10; 
end

varNames = {'i','t','x','y','v_x','v_y'};
Runge_Kutta=table((0:1:dis-1)',t(1:dis),x(1:dis),y(1:dis),vx(1:dis),vy(1:dis),'VariableNames',varNames);
disp(Runge_Kutta)

figure 
plot(x,y,'*')
grid('on')
ylim([-2*b 1.5*b])
hold('on')
C = [-c 0] ;                 
th = linspace(0,2*pi) ; 
xe = C(1)+a*cos(th) ; 
ye = C(2)+b*sin(th) ; 
plot(xe,ye,'r')
axis equal
plot(0,0,'o','Color','r')
legend('3nd order Runge-Kutta method','Exact solution','Focus S','Location','south')


disp(' ')
fprintf('   ======================= \n')
fprintf('   ==== Error analysis === \n')
fprintf('   ======================= \n \n')

res=[t w_RK];                       % Results of the 3rd order Runge-Kutta solver
u=NaN(n+1,1);                       % Eccentric anomalies 
for i=1:n+1
   Ma=2*pi/T*t(i);                  % Mean anomaly 
   diff=1;                          % Resolution of the Kepler's equation 
   u_0=Ma/(1-e);                    % Initialisation u0
   while abs(diff)>0.0001           % Precision criteria for the resolution 
       u_1=e*sin(u_0)+Ma;           % Recurrence equation
       diff=u_1-u_0;                % Difference between u_1 and u_0
       u_0=u_1;                     % Previous value
   end                              % End of the resolution of the Kepler's equation 
   u(i)=u_1;                        % Excentric anomaly
end

v=atan(sqrt((1+e)/(1-e))*tan(u/2))*2;  % True anomaly
x_true=a*(cos(u)-e);                   % True position x
y_true=a*sin(u)*sqrt(1-e^2);           % True position y
re_x=abs(x-x_true)./abs(x_true)*100;   % Relative error on x
re_y=abs(y-y_true)./abs(y_true)*100;   % Relative error on y

varNames = {'i','t','x','y','x true','y true','Relative error x (%)','Relative error y (%)'};
Error_analysis=table((0:1:dis-1)',t(1:dis),x(1:dis),y(1:dis),x_true(1:dis),y_true(1:dis),re_x(1:dis),re_y(1:dis),'VariableNames',varNames);
disp(Error_analysis)

fprintf('   Step size = %1.4f \n',h)
fprintf('   Mean of relative errors on x (percentage) = %1.4f \n',nanmean(re_x))
fprintf('   Mean of relative errors on y (percentage) = %1.4f \n',nanmean(re_y))
fprintf('   Mean of absolute errors on x = %1.6f \n',nanmean(abs(x-x_true)))
fprintf('   Mean of absolute errors on y = %1.6f \n',nanmean(abs(y-y_true)))

