clear all
clc

earthmain

function earthmain 
mu = 1.32712440018e20; % m^3/s^2
r = 149597870700; % m  (= 1 AU)
vcirc = sqrt(mu/r); % m/s (circular velocity for the given value of r)
y0 = [r; 0; 0; 0; vcirc; 0]; % m and m/s (start on x-axis with +y velocity direction)
tspan=[0, 365.25*86400]; % One year time span
options = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,y]=ode45(@(t,y)earthtwo(t,y,mu),tspan,y0,options);
% for i = 1:length (t) 
%     fprintf ('%2i %7.5f %7.5f \n',i,t(i),y(i)); 
% end
figure;
plot(y(:,1),y(:,2),'*'); % Plot only the x-y plane since that is how we set things up
axis square
grid on
title('The solution will give us the orbit of a body in two dimensions.'); 
xlabel('x (m)'); 
ylabel('y (m)');
end
function f= earthtwo(t,y,mu) 
rx = y(1); 
ry = y(2); 
rz = y(3); 
vx = y(4); 
vy = y(5); 
vz = y(6); 
r = sqrt(rx.^2+ry.^2+rz.^2);
f=zeros(6,1); 
f(1) = vx; 
f(2) = vy; 
f(3) = vz; 
f(4) = -mu*rx/r.^3; 
f(5) = -mu*ry/r.^3; 
f(6) = -mu*rz/r.^3; 
end
