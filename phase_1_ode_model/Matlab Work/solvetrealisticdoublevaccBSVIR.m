%Solves the realistic double vaccination schedual of the BSVIR model using
%the fourth stage Runge-Kutta method. Storing the partially age stratified
%model as a matrix. It then plots the number of susceptible and infectious
%individuals.
%n - the number of years to run the system
%lowvage - the age at which first vaccination schedual is applied
%highvage - the age at which second vaccination schedual is applied
%vyear - year in which vaccination is introduced.
%p - Proportion of population vaccinated
%u0 - intial number individuals in each group. u0=[br0;bv0;s0;v0;i0;r0]
%population - the matrix which stores age information and status
%information.
%t - current time
%h - time step
%c - counter
%introvaccstep - the step at which vaccine is introduced.
%lowvagegroup - the low age group at which vaccination is applied. 
%highvagegroup - the high age group at which vaccination is applied. 
%l - length of t
%U - tracks the total number in each group over time.
%***Due to computation required to calulate, it is effective to store 
%values and reference them when required.***
%lossofunvaccimmunityvalues - Store proportion of infants of unvaccinated
%mothers that lose immunity between time steps.
%lossofvaccimmunityvalues - Store proportion of infants of vaccinated
%mothers that lose immunity between time steps.
function y = solvetrealisticdoublevaccBSVIR(u0,p,lowvage,highvage,vyear,n)
%Intial Values
%Br0 = 0;
%Bv0 = 0;
%s0 = 6400000;
%v0 = 0;
%i0 = 600000;
%r0 = 53000000;
%u0 = [0;0;s0;0;i0;r0];
c = 1;
h = 1/365;
t = [0:h:n];
l = length(t);
introvaccstep = vyear/h;
lowvagegroup = lowvage/h;
highvagegroup = highvage/h;
%Population matrix stores information about infants until they reach
%second vaccination age. After that they are added to the general group,
%the last column of the matrix.
population = zeros(6,highvagegroup+1);
population(:,end) = u0;
lossofunvaccimmunityvalues = zeros (1,highvagegroup);
for i=1:length(lossofunvaccimmunityvalues)
   lossofunvaccimmunityvalues(i) =  lossofunvaccinatedimmunity(i*h,(i+1)*h);
end
lossofvaccimmunityvalues = zeros (1,highvagegroup);
for i=1:length(lossofvaccimmunityvalues)
   lossofvaccimmunityvalues(i) =  lossofvaccinatedimmunity(i*h,(i+1)*h); 
end
U = zeros(6,l);
while ((c<= l)&&(c<=introvaccstep))
    %Runge-Kutta Method for new birth group
    kbirth1 = birthBSVIR(t(c),population(:,end));
    kbirth2 = birthBSVIR(t(c)+h/2,population(:,end)+(h*kbirth1/2));
    kbirth3 = birthBSVIR(t(c)+h/2,population(:,end)+(h*kbirth2/2));
    kbirth4 = birthBSVIR(t(c)+h,population(:,end)+(h*kbirth3));
    births = (h*(kbirth1+(2*kbirth2)+(2*kbirth3)+kbirth4)/6);
    %Runge-Kutta Method for change in population dynamics
    k1 = matrixBSVIR(t(c),population);
    k2 = matrixBSVIR(t(c)+h/2,population+(h*k1/2));
    k3 = matrixBSVIR(t(c)+h/2,population+(h*k2/2));
    k4 = matrixBSVIR(t(c)+h,population+(h*k3));
    population = population + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    for i=1:(length(population)-1) 
        brnum = population(1,i);
        bvnum = population(2,i);
        population(1,i) = lossofunvaccimmunityvalues(i)*brnum;
        population(2,i) = lossofvaccimmunityvalues(i)*bvnum;
        population(3,i) = population(3,i)+((1-lossofunvaccimmunityvalues(i))*brnum)+((1-lossofvaccimmunityvalues(i))*bvnum);
    end
    %Shift data
    population = [births,population(:,1:(end-2)),(population(:,end-1)+population(:,end))];
    U(:,c) = sum(population,2);
    c = c+1;
end
while (c<=l)
    %Runge-Kutta Method for new birth group
    kbirth1 = birthBSVIR(t(c),population(:,end));
    kbirth2 = birthBSVIR(t(c)+h/2,population(:,end)+(h*kbirth1/2));
    kbirth3 = birthBSVIR(t(c)+h/2,population(:,end)+(h*kbirth2/2));
    kbirth4 = birthBSVIR(t(c)+h,population(:,end)+(h*kbirth3));
    births = (h*(kbirth1+(2*kbirth2)+(2*kbirth3)+kbirth4)/6);
    %Apply Vaccination
    %Lower age group
    population(4,lowvagegroup) = p*population(3,lowvagegroup);
    population(3,lowvagegroup) = (1-p)*population(3,lowvagegroup);
    %Higher age group
    population(4,highvagegroup) = p*population(3,highvagegroup);
    population(3,highvagegroup) = (1-p)*population(3,highvagegroup);
    %Runge-Kutta Method for change in population dynamics
    k1 = matrixBSVIR(t(c),population);
    k2 = matrixBSVIR(t(c)+h/2,population+(h*k1/2));
    k3 = matrixBSVIR(t(c)+h/2,population+(h*k2/2));
    k4 = matrixBSVIR(t(c)+h,population+(h*k3));
    population = population + (h*(k1+(2*k2)+(2*k3)+k4)/6);
    for i=1:(length(population)-1) 
        brnum = population(1,i);
        bvnum = population(2,i);
        population(1,i) = lossofunvaccimmunityvalues(i)*brnum;
        population(2,i) = lossofvaccimmunityvalues(i)*bvnum;
        population(3,i) = population(3,i)+((1-lossofunvaccimmunityvalues(i))*brnum)+((1-lossofvaccimmunityvalues(i))*bvnum);
    end
    %Shift data
    population = [births,population(:,1:(end-2)),(population(:,end-1)+population(:,end))];
    U(:,c) = sum(population,2);
    c = c+1;
end
y = U;
%Display the resulting susceptible and infectious individuals.
%S = U(3,:);
%I = U(5,:);
%plot(t,S,'r',t,I,'b')
%legend('Susceptable','Infectious')
%x_label = xlabel('Years');
%set(x_label,'FontSize',14);
%y_label = ylabel('Current number of individuals');
%set(y_label,'FontSize',14);