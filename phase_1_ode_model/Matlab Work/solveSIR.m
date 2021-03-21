%Solves the simple SIR model using the fourth stage Runge-Kutta method.
%It then plots the number of susceptible and infectious individuals.
%n - the number of years to run the system
%s0 - intial number susceptable in population
%i0 - intial number infected in population
%r0 - intial number recovered in population
%t - current time (vector)
%h - time step
%c - counter
%l - length of t
%U - Stores the number of individuals in each group at each time step.
function y = solveSIR(n)
%Intial Values
s0 = 6400000;
i0 = 600000;
r0 = 53000000;
u0 = [s0;i0;r0];
u = u0;
c = 1;
h = 1/365;
t = [0:h:n];
l = length(t);
U = zeros(3,l);
while c<= l
    %Runge-Kutta Method
    k1 = SIR(t(c),u);
    k2 = SIR(t(c)+h/2,u+(h*k1/2));
    k3 = SIR(t(c)+h/2,u+(h*k2/2));
    k4 = SIR(t(c)+h,u+(h*k3));
    u = u + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    U(:,c) = u;
    c = c+1;
end
y = U;
%Display the resulting susceptible and infectious individuals.
S = U(1,:);
I = U(2,:);
plot(t,S,'r',t,I,'b')
h_legend=legend('Susceptible','Infectious');
set(h_legend,'FontSize',14);
x_label = xlabel('Years');
set(x_label,'FontSize',18);
y_label = ylabel('Current number of individuals');
set(y_label,'FontSize',16);