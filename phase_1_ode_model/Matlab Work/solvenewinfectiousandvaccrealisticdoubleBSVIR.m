%Solves the realistic (double vaccination) BSVIR model using the fourth
%stage Runge-Kutta method. Storing the partially age stratified model 
%as a matri. Returning the numbers of new infectious and vaccinated
%individuals at each time step.
%n - the number of years to run the system
%lowvage - the first age at which vaccination schedual is applied
%highvage - the second age at which vaccination schedual is applied
%vyear - year in which vaccination is introduced.
%p - Proportion of population vaccinated
%population - the matrix which input of the population.
%information.
%t - current time
%h - time step
%c - counter
%introvaccstep - the step at which vaccine is introduced.
%lowvagegroup - the low age group at which vaccination is applied. 
%highvagegroup - the high age group at which vaccination is applied.
%l - length of t
%lpop - length of the population matrix
%U - tracks the total number in each group over time.
%newinfectiousindivuals - stores the number of newly infected individuals
%as time progresses.
%newvaccinatedindivuals - stores the number of newly vaccinated individuals
%as time progresses.
%***Due to computation required to calulate, it is effective to store 
%values and reference them when required.***
%lossofunvaccimmunityvalues - Store proportion of infants of unvaccinated
%mothers that lose immunity between time steps.
%lossofvaccimmunityvalues - Store proportion of infants of vaccinated
%mothers that lose immunity between time steps.
function y = solvenewinfectiousandvaccrealisticdoubleBSVIR(population,p,lowvage,highvage,vyear,n)
c = 1;
h = 1/365;
t = [0:h:n];
l = length(t);
introvaccstep = vyear/h;
lowvagegroup = lowvage/h;
highvagegroup = highvage/h;
%Population matrix stores information about infants until they reach
%vaccination age. After that they are added to the general group, the last
%column of the matrix.
lpop = length(population);
lossofunvaccimmunityvalues = zeros (1,lpop-1);
for i=1:length(lossofunvaccimmunityvalues)
   lossofunvaccimmunityvalues(i) =  lossofunvaccinatedimmunity(i*h,(i+1)*h);
end
lossofvaccimmunityvalues = zeros (1,lpop-1);
for i=1:length(lossofvaccimmunityvalues)
   lossofvaccimmunityvalues(i) =  lossofvaccinatedimmunity(i*h,(i+1)*h); 
end
U = zeros(6,l);
newinfectiousindividuals = zeros(1,l);
newvaccinatedindividuals = zeros(1,l);
while ((c<= l)&&(c<=introvaccstep))
    %Runge-Kutta Method for new birth group
    kbirth1 = birthBSVIR(t(c),population(:,end));
    kbirth2 = birthBSVIR(t(c)+h/2,population(:,end)+(h*kbirth1/2));
    kbirth3 = birthBSVIR(t(c)+h/2,population(:,end)+(h*kbirth2/2));
    kbirth4 = birthBSVIR(t(c)+h,population(:,end)+(h*kbirth3));
    births = (h*(kbirth1+(2*kbirth2)+(2*kbirth3)+kbirth4)/6);
    %Runge-Kutta Method for the number of new infectious individuals
    k1 = newinfectiousrealisticBSVIR(t(c),sum(population,2));
    k2 = newinfectiousrealisticBSVIR(t(c)+h/2,sum(population,2)+(h*k1/2));
    k3 = newinfectiousrealisticBSVIR(t(c)+h/2,sum(population,2)+(h*k2/2));
    k4 = newinfectiousrealisticBSVIR(t(c)+h,sum(population,2)+(h*k3));
    delta = (h*(k1+(2*k2)+(2*k3)+k4)/6);
    newinfectiousindividuals(c) = delta (5);
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
    %Runge-Kutta Method for the number of new infectious individuals
    k1 = newinfectiousrealisticBSVIR(t(c),sum(population,2));
    k2 = newinfectiousrealisticBSVIR(t(c)+h/2,sum(population,2)+(h*k1/2));
    k3 = newinfectiousrealisticBSVIR(t(c)+h/2,sum(population,2)+(h*k2/2));
    k4 = newinfectiousrealisticBSVIR(t(c)+h,sum(population,2)+(h*k3));
    delta = (h*(k1+(2*k2)+(2*k3)+k4)/6);
    newinfectiousindividuals(c) = delta (5);
    %Apply Vaccination
    population(4,lowvagegroup) = p*population(3,lowvagegroup);
    population(3,lowvagegroup) = (1-p)*population(3,lowvagegroup);
    population(4,highvagegroup) = p*population(3,highvagegroup)+population(4,highvagegroup);
    population(3,highvagegroup) = (1-p)*population(3,highvagegroup);
    %Add the number of vaccinated individuals.
    newvaccinatedindividuals(c)= p*(sum(population(:,lowvagegroup))+sum(population(:,highvagegroup)));
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
y = [newinfectiousindividuals;newvaccinatedindividuals];
%Display the number new infectious and vaccinated individuals at each time
%step.
%plot(t,newvaccinatedindividuals,'r',t,newinfectiousindividuals,'b')
%legend('Vaccinated','Infectious')
%x_label = xlabel('Years');
%set(x_label,'FontSize',14);
%y_label = ylabel('Number of new individuals');
%set(y_label,'FontSize',14);