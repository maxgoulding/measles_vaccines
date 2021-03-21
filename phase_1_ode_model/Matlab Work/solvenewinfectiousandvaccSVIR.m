%Solves the new number of infectious and vaccinated individuals in 
%the simple SVIR model using the fourth-stage Runge-Kutta method. It
%then solves the simple SVIR model using the fourth stage Runge-Kutta
%method to track overall system behaviour. It then plots the new number
%infectious and vaccinated individuals.
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
%introvaccstep - the step at which vaccination is introduced
%l - length of t
%U - Stores the number of individuals in each group at each time step.
%newinfectiousindivuals - stores the number of newly infected individuals
%as time progresses.
%newvaccinatedindivuals - stores the number of newly vaccinated individuals
%as time progresses.
%delta - temporarily stores the changes in the newinfectiousandvaccSVIR
%differential model. i.e. the newly infected and vaccinated individuals.
function y = solvenewinfectiousandvaccSVIR(u0,p,vyear,n)
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
newinfectiousindividuals = zeros(1,l);
newvaccinatedindividuals = zeros(1,l);
while ((c<= l)&&(c<=introvaccstep))
    %Runge-Kutta Method for the number of new infectious individuals
    k1 = newinfectiousandvaccSVIR(0,t(c),u);
    k2 = newinfectiousandvaccSVIR(0,t(c)+h/2,u+(h*k1/2));
    k3 = newinfectiousandvaccSVIR(0,t(c)+h/2,u+(h*k2/2));
    k4 = newinfectiousandvaccSVIR(0,t(c)+h,u+(h*k3));
    delta = (h*(k1+(2*k2)+(2*k3)+k4)/6);
    newvaccinatedindividuals(c) = delta(2);
    newinfectiousindividuals(c) = delta(3);
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
    %Runge-Kutta Method for the number of new infectious individuals
    k1 = newinfectiousandvaccSVIR(p,t(c),u);
    k2 = newinfectiousandvaccSVIR(p,t(c)+h/2,u+(h*k1/2));
    k3 = newinfectiousandvaccSVIR(p,t(c)+h/2,u+(h*k2/2));
    k4 = newinfectiousandvaccSVIR(p,t(c)+h,u+(h*k3));
    delta = (h*(k1+(2*k2)+(2*k3)+k4)/6);
    newvaccinatedindividuals(c) = delta(2);
    newinfectiousindividuals(c) = delta(3);
    %Runge-Kutta Method
    k1 = SVIR(p,t(c),u);
    k2 = SVIR(p,t(c)+h/2,u+(h*k1/2));
    k3 = SVIR(p,t(c)+h/2,u+(h*k2/2));
    k4 = SVIR(p,t(c)+h,u+(h*k3));
    u = u + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    U(:,c) = u;
    c = c+1;
end
y = [newinfectiousindividuals;newvaccinatedindividuals];
%Display the number new infectious and vaccinated individuals at each time
%step.
plot(t,newvaccinatedindividuals,'r',t,newinfectiousindividuals,'b')
h_legend=legend('Vaccinated','Infectious');
set(h_legend,'FontSize',14);
x_label = xlabel('Years');
set(x_label,'FontSize',18);
y_label = ylabel('Number of new individuals');
set(y_label,'FontSize',16);