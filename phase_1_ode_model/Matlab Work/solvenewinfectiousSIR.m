%Solves the new number of infectious individuals in the simple SIR model
%using the fourth-stage Runge-Kutta method. It then solves the simple SIR
%model using the fourth stage Runge-Kutta method to track overall system
%behaviour. It then plots the new number infectious individuals.
%n - the number of years to run the system
%s0 - intial number susceptable in population
%i0 - intial number infected in population
%r0 - intial number recovered in population
%t - current time (vector)
%h - time step
%c - counter
%l - length of t
%U - Stores the number of individuals in each group at each time step.
%newinfectiousindivuals - stores the number of newly infected individuals
%as time progresses.
%delta - temporarily stores the changes in the newinfectiousSIR
%differential model. i.e. the newly infected individuals.
function y = solvenewinfectiousSIR(n)
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
newinfectiousindividuals = zeros(1,l);
while c<= l
    %Runge-Kutta Method for the number of new infectious individuals
    k1 = newinfectiousSIR(t(c),u);
    k2 = newinfectiousSIR(t(c)+h/2,u+(h*k1/2));
    k3 = newinfectiousSIR(t(c)+h/2,u+(h*k2/2));
    k4 = newinfectiousSIR(t(c)+h,u+(h*k3));
    delta = (h*(k1+(2*k2)+(2*k3)+k4)/6);
    newinfectiousindividuals(c) = delta(2);
    %Runge-Kutta Method for the SIR model
    k1 = SIR(t(c),u);
    k2 = SIR(t(c)+h/2,u+(h*k1/2));
    k3 = SIR(t(c)+h/2,u+(h*k2/2));
    k4 = SIR(t(c)+h,u+(h*k3));
    u = u + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    U(:,c) = u;
    c = c+1;
end
y = newinfectiousindividuals;
%Display the resulting new infectious individuals.
plot(t,newinfectiousindividuals,'b')
x_label = xlabel('Years');
set(x_label,'FontSize',18);
y_label = ylabel('Number of new individuals');
set(y_label,'FontSize',16);