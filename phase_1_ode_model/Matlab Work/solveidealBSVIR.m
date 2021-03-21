%Solves the ideal BSVIR model using the fourth stage Runge-Kutta method.
%It then plots the number of susceptible and infectious individuals.
%u0 - intial distribution of population
%p - Proportion of population vaccinated
%vyear - year in which vaccination is introduced
%n - the number of years to run the system
%t - current time
%h - time step
%c - counter
%introvaccstep - the step at which vaccine is introduced
%l - length of t
%U - Store of population groups at each time step
function y = solveidealBSVIR(u0,p,vyear,n)
u = u0;
c = 1;
h = 1/365;
t = [0:h:n];
introvaccstep = vyear/h;
l = length(t);
U = zeros(6,l);
while ((c<= l)&&(c<=introvaccstep))
    %Runge-Kutta Method
    k1 = idealBSVIR(0,t(c),u);
    k2 = idealBSVIR(0,t(c)+h/2,u+(h*k1/2));
    k3 = idealBSVIR(0,t(c)+h/2,u+(h*k2/2));
    k4 = idealBSVIR(0,t(c)+h,u+(h*k3));
    u = u + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    U(:,c)=u;
    c = c+1;
end
while (c<=l)
    %Runge-Kutta Method
    k1 = idealBSVIR(p,t(c),u);
    k2 = idealBSVIR(p,t(c)+h/2,u+(h*k1/2));
    k3 = idealBSVIR(p,t(c)+h/2,u+(h*k2/2));
    k4 = idealBSVIR(p,t(c)+h,u+(h*k3));
    u = u + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    U(:,c) = u;
    c = c+1;
end
y = U;
%Display the resulting susceptible and infectious individuals.
%S = U(3,:);
%I = U(5,:);
%plot(t,S,'r',t,I,'b')
%h_legend=legend('Susceptible','Infectious');
%set(h_legend,'FontSize',14);
%x_label = xlabel('Years');
%set(x_label,'FontSize',18);
%y_label = ylabel('Current number of individuals');
%set(y_label,'FontSize',16);