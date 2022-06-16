function Chemotaxis_Presentation_Data_Generation

% Specify alpha values at which to run the code at. To run the base
% simulation, set 'a = 1'
for a = [0.125:0.125:2]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inital Conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N      = 400; % Number of parts interval [0,1] is broken into
h      = 1/N; % Length of eacch interval
ep1    = 0.01; % Coeffificient of diffusion of v
ep2    = ep1; % Coefficient of advection of v
coeff1 = 1;   % Logistic term coefficient
T      = 2;   % Total system runtime
alpha  = a;   % Logistic growth exponent

write_to_csv = 1000; % Write the output to a .csv file once every __ outputs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% u_t + (uv)_x = u_xx + coeff1 * u(1-u)
% v_t + u_x = ep1 * v_xx + ep2 * (v^2)_x

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filenames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Instructions for how to format filenames, %f for floating point numbers.
% Create a file for both u and v that shows inputs
formatSpec_u = 'u_%.2f_%.2f.csv';
formatSpec_v = 'v_%.2f_%.2f.csv';

% Create the filename strings with above instructions
filename_u = sprintf(formatSpec_u,T,alpha);
filename_v = sprintf(formatSpec_v,T,alpha);

% Create a filename for timescale results
formatSpecResults = 'results_%.2f_%.2f.csv';
filename_results = sprintf(formatSpecResults,T,alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding Initial Conditions to Results .csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writematrix(sprintf('%.0f',N),filename_results,'WriteMode','append')
writematrix(sprintf('%f',h),filename_results,'WriteMode','append')
writematrix(sprintf('%f',ep1),filename_results,'WriteMode','append')
writematrix(sprintf('%f',ep2),filename_results,'WriteMode','append')
writematrix(sprintf('%f',coeff1),filename_results,'WriteMode','append')
writematrix(sprintf('%f',alpha),filename_results,'WriteMode','append')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Conditions and Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Stop the calculation if == 1
% Used for early termination of the program
stopv=0;
stopu=0;

% Create the initial x-values for u and v
% Both offcenter and different dimensions
% See Dr. Fuster's notes
xv = h*(0:1:N)';
xu = h*( 1/2 : 1 : N-1/2)';

% Initial condition functions for u and v
u0 = 5/6 + erf((xu-1/2)*8) / pi;
v0 = 10 * exp(-(xv-1/2).^2*18) / 6;

% Correct v0 so that v0(0) and v0(1) == 0
v0 = v0-v0(1);

% Calculate the average of the function over the interval
ave = sum(u0*h);

% Create an array of ones
one=ones(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write initial data to a .csv file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

writematrix(reshape(u0, [1,N]),filename_u)
writematrix(reshape(v0, [1,N+1]),filename_v)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Precalculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inv_h  = 1 / h;
inv_h2 = 1 / (h ^ 2);
inv_2h = 1 / (2 * h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Applying initial values to system and final parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = u0;
v = v0;
dt = h^2 * 0.5; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Running the simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp('epsilon='), disp(ep1), disp('logistic r='), disp(coeff1)
%plot(xu,u0,xv,v0,'LineWidth',2.0),axis([0 1 0 2]),legend('u','v'), xlabel('x'), pause(0.05)  %plots initial functions, axis is [0,1] (we calculate from 0 to 1.1??why?, time step is 0.5 

% Finds how many timesteps fit in the runtime and makes sure it is an integer.
for k = 1:floor(T/dt)

% Establish Neumann boundary conditions of u_x = 0 for u
uAUG = [u(1);u;u(end)];
%uAUG(2:N+2)
uv  = (uAUG(2:N+2)+uAUG(1:N+1)).* v * 0.5;  %calculates all needed/ derivatives (2nd order precision)
uvx = (uv(2:N+1)-uv(1:N))* inv_h;
uxx = (uAUG(3:N+2)-2*u+uAUG(1:N))*inv_h2;

vx=( v(3:N+1) - v(1:N-1) )*inv_2h;
v2 = v.^2;
v2x = ( v2(3:N+1) - v2(1:N-1) )*inv_2h;
ux  = ( uAUG(3:N+1) - uAUG(2:N) )*inv_h;
vxx = (v(1:N-1)-2*v(2:N)+v(3:N+1))*inv_h2;

unp1 = u + dt*( uxx - uvx + coeff1*u.*((1-u).^a) ); %u_t= (unp1- u)/ dt  %solves the diff equation
vnp1 = v(2:N) + dt*( ep1*vxx + ep2*v2x - ux ); %v_t= (vnp1 - v)/dt

%new variables to compare/store
%vxp1= (vnp1(3:N-1)-vnp1(1:N-3) )/(2*h);
uAUGp1=[unp1(1);unp1;unp1(end)];
%uxp1=( uAUGp1(3:N+1) - uAUGp1(2:N) )/(h);
uxxp1= (uAUGp1(3:N+2)-2*unp1+uAUGp1(1:N))*inv_h2;
vxxp1=(vnp1(1:N-3)-2*vnp1(2:N-2)+vnp1(3:N-1))*inv_h2;

t=k*dt;
%sum(abs(vxp1).^2)/(N-1)
%sum(abs(vx(2:N-1)).^2)/(N-2)
%sqrt(sum(abs(uxp1).^2))/N
%sqrt(sum(abs(vxp1).^2))/N
if abs((sum(abs(unp1-ave).^2))/(N))<0.01 && abs((sum(abs(u-ave).^2))/(N))>0.01
    disp('initial time close to average')
    writematrix(sprintf('initial time close to average = %f',t),filename_results,'WriteMode','append')
    tinitial=t
    %figure(1);  plot(xu,u,xv,v), xlabel('x'),legend('u','v'), axis([0 1 -.1 1.1]),title(['t=',num2str(k*dt),'r=',num2str(coeff1), 'eps =',num2str(ep1)]), pause(2) 
end
if abs((sum(abs(unp1-ave).^2))/(N))>0.01 && abs((sum(abs(u-ave).^2))/(N))<0.01

     disp('final time close to average')
    writematrix(sprintf('final time close to average = %f',t),filename_results,'WriteMode','append')
    tfinal=t;
    %figure(2); plot(xu,u,xv,v), xlabel('x'),legend('u','v'), axis([0 1 -.1 1.1]),title(['t=',num2str(k*dt),'r=',num2str(coeff1), 'eps =',num2str(ep1)]), pause(2) 
end
if abs((sum(abs(unp1-one).^2))/(N))<0.01 && abs((sum(abs(u-one).^2))/(N))>=0.01

    disp('initial time close to carrying capacity')
    writematrix(sprintf('initial time close to carrying capacity = %f',t),filename_results,'WriteMode','append')
    tcarrying=t
    %figure(3); plot(xu,u,xv,v), xlabel('x'),legend('u','v'), axis([0 1 -.1 1.1]),title(['t=',num2str(k*dt),'r=',num2str(coeff1), 'eps =',num2str(ep1)]), pause(2) 
end
if  sum(abs(ux).^2)/(N-1)<0.01 && stopu==0

    disp('relaxation for u time')
    writematrix(sprintf('relaxation for u time = %f',t),filename_results,'WriteMode','append')
    trelaxu=t
    stopu=1;
    %figure(4); plot(xu,u,xv,v), xlabel('x'),legend('u','v'), axis([0 1 -.1 1.1]),title(['t=',num2str(k*dt),'r=',num2str(coeff1), 'eps =',num2str(ep1)]), pause(2) 
    
end
if sum(abs(vx).^2)/(N-1)<0.01 && stopv==0

    disp('relaxation for v time')
    writematrix(sprintf('relaxation for v time = %f',t),filename_results,'WriteMode','append')
    trelaxv=t
    stopv=1;
     %figure(5); plot(xu,u,xv,v), xlabel('x'),legend('u','v'), axis([0 1 -.1 1.1]),title(['t=',num2str(k*dt),'r=',num2str(coeff1), 'eps =',num2str(ep1)]), pause(0.1) 
end
if sum(abs(vxxp1).^2)/N<0.01 && sum(abs(vxx.^2))/N>0.01 

    disp('relaxation for vxx time')
    writematrix(sprintf('relaxation for vxx time = %f',t),filename_results,'WriteMode','append')
    trelaxvxx=t
    
     %figure(5); plot(xu,u,xv,v), xlabel('x'),legend('u','v'), axis([0 1 -.1 1.1]),title(['t=',num2str(k*dt),'r=',num2str(coeff1), 'eps =',num2str(ep1)]), pause(0.1) 
end
u    = unp1; % gets the new value calculated for what u is  
v    = [0;vnp1;0]; % gets the new value calculated for what v is but keeps it with dirichlet b.c. at each step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the output to the .csv file once every N iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(k,write_to_csv) == 1    
    writematrix(reshape(real(u), [1,N]),filename_u,'WriteMode','append')
    writematrix(reshape(real(v), [1,N+1]),filename_v,'WriteMode','append')
    %plot(xu,u,xv,v,'LineWidth',2.0), xlabel('x'),legend('u','v'), axis([0 1 -.1 2]),title(['t=',num2str(k*dt),'r=',num2str(coeff1), 'eps =',num2str(ep1)]), pause(0.05) 
  
end

end

end
