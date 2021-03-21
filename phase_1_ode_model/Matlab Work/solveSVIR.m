%Solves the SVIR model using the fourth stage Runge-Kutta method.
%It then plots the number of susceptible and infectious individuals.
%n - the number of years to run the system
%vyear - year in which vaccination is introduced.
%p - Proportion of newborns vaccinate
%s0 - intial number susceptable in population
%v0 - intial number vaccinated in population
%i0 - intial number infected in population
%r0 - intial number recovered in population
%t - current time
%h - time step
%c - counter
%introvaccstep - the step at which vaccine is introduced.
%l - length of t
function y = solveSVIR(u0,p,vyear,n)
%Intial Values
%s0 = 6400000;
%v0 = 0;
%i0 = 600000;
%r0 = 53000000;
%u0 = [s0;v0;i0;r0];
u = u0;
c = 1;
h = 1/365;
t = [0:h:n];
introvaccstep = vyear/h;
l = length(t);
U = zeros(4,l);
while ((c<= l)&&(c<=introvaccstep))
    %Runge-Kutta Method
    k1 = SVIR(0,t(c),u);
    k2 = SVIR(0,t(c)+h/2,u+(h*k1/2));
    k3 = SVIR(0,t(c)+h/2,u+(h*k2/2));
    k4 = SVIR(0,t(c)+h,u+(h*k3));
    u = u + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    U(:,c) = u;
    c = c+1;
end
while (c<=l)
    %Runge-Kutta Method
    k1 = SVIR(p,t(c),u);
    k2 = SVIR(p,t(c)+h/2,u+(h*k1/2));
    k3 = SVIR(p,t(c)+h/2,u+(h*k2/2));
    k4 = SVIR(p,t(c)+h,u+(h*k3));
    u = u + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    U(:,c) = u;
    c = c+1;
end
y = U;
%Display the resulting susceptible and infectious individuals.
S = U(1,:);
I = U(3,:);
%plot(t,S,'r',t,I,'b')
%h_legend=legend('Susceptible','Infectious');
%set(h_legend,'FontSize',14);
%x_label = xlabel('Years');
%set(x_label,'FontSize',18);
%y_label = ylabel('Current number of individuals');
%set(y_label,'FontSize',16);