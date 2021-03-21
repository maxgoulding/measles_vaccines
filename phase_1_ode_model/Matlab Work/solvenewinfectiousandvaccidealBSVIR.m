%Solves the new number of infectious and vaccinated individuals in 
%the ideal BSVIR model using the fourth-stage Runge-Kutta method. It
%then solves the ideal BSVIR model using the fourth stage Runge-Kutta
%method to track overall system behaviour. It then plots the new number
%infectious and vaccinated individuals.
%u0 - initial values
%p - Proportion of population vaccinated
%vyear - year in which vaccination is introduced.
%n - the number of years to run the system
%t - current time
%h - time step
%c - counter
%introvaccstep - the step at which vaccine is introduced
%l - length of t
%U - Stores the number of individuals in each group at each time step.
%newinfectiousindivuals - stores the number of newly infected individuals
%as time progresses.
%newvaccinatedindivuals - stores the number of newly vaccinated individuals
%as time progresses.
%delta - temporarily stores the changes in the newinfectiousandvaccSVIR
%differential model. i.e. the newly infected and vaccinated individuals.
function y = solvenewinfectiousandvaccidealBSVIR(u0,p,vyear,n)
u = u0;
c = 1;
h = 1/365;
t = [0:h:n];
introvaccstep = vyear/h;
l = length(t);
U = zeros(6,l);
newinfectiousindividuals = zeros(1,l);
newvaccinatedindividuals = zeros(1,l);
while ((c<= l)&&(c<=introvaccstep))
    %Runge-Kutta Method for the number of new infectious individuals
    k1 = newinfectiousandvaccidealBSVIR(0,t(c),u);
    k2 = newinfectiousandvaccidealBSVIR(0,t(c)+h/2,u+(h*k1/2));
    k3 = newinfectiousandvaccidealBSVIR(0,t(c)+h/2,u+(h*k2/2));
    k4 = newinfectiousandvaccidealBSVIR(0,t(c)+h,u+(h*k3));
    delta = (h*(k1+(2*k2)+(2*k3)+k4)/6);
    newvaccinatedindividuals(c) = delta(4);
    newinfectiousindividuals(c) = delta(5);
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
    %Runge-Kutta Method for the number of new infectious and vaccinated
    %individuals
    k1 = newinfectiousandvaccidealBSVIR(p,t(c),u);
    k2 = newinfectiousandvaccidealBSVIR(p,t(c)+h/2,u+(h*k1/2));
    k3 = newinfectiousandvaccidealBSVIR(p,t(c)+h/2,u+(h*k2/2));
    k4 = newinfectiousandvaccidealBSVIR(p,t(c)+h,u+(h*k3));
    delta = (h*(k1+(2*k2)+(2*k3)+k4)/6);
    newvaccinatedindividuals(c) = delta(4);
    newinfectiousindividuals(c) = delta(5);
    %Runge-Kutta Method
    k1 = idealBSVIR(p,t(c),u);
    k2 = idealBSVIR(p,t(c)+h/2,u+(h*k1/2));
    k3 = idealBSVIR(p,t(c)+h/2,u+(h*k2/2));
    k4 = idealBSVIR(p,t(c)+h,u+(h*k3));
    u = u + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    U(:,c) = u;
    c = c+1;
end
y = [newinfectiousindividuals;newvaccinatedindividuals];
%Display the number new infectious and vaccinated individuals at each time
%step.
%plot(t,newvaccinatedindividuals,'r',t,newinfectiousindividuals,'b')
%legend('Vaccinated','Infectious')
%x_label = xlabel('Years');
%set(x_label,'FontSize',14);
%y_label = ylabel('Number of new individuals');
%set(y_label,'FontSize',14);