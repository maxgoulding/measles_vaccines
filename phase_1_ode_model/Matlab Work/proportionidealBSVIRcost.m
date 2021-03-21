%Determines the cost of vaccinating varying proportions of the population
%against measles, using the ideal BSVIR model as a basis for the number of
%infections and vaccinations occured. The system is run for 10 years
%without vaccination, vaccination is then applied at 85% for 20 years to
%give a mass-vaccinated population. The system is then run for an
%additional 50 years to determine the long term behaviour. The number of
%new infectious individuals and vaccinations that occured is then summed to
%determine those incidence which incur cost. This is then multiplied by the
%associated cost to give total cost of the 50 year period from mass-vaccination.
%Parameters
%costpermeaslescase - the cost of a measles case in US $ (2001) to the UK.
%costpervaccination - the cost of measles vaccination in US $ (2001) to the UK.
%associatedcostpervaccination - the additional associated cost of each
%vaccination, e.g. any side effects of the vaccine.
%infections - the number of new infections that occur at each proportion of
%the population vaccinated.
%vaccinations - the number of vaccinations that occur at each proportion of
%the population vaccinated.
%intromassvaccyear - the year in which mass vaccination is introduced.
%massvaccperiod - the time which mass vaccination occured for. Years
%massvaccprop - the proportion of the population that is vaccinated during
%the mass period.
%s0 - intial number susceptable in population
%v0 - intial number vaccinated in population
%i0 - intial number infected in population
%r0 - intial number recovered in population
%U - store of the values for each group obtained during the mass
%vaccination period.
%postmassvacc - Is a vector giving the 'intial' values after the mass
%vaccination period. This is the intial values for the next 50 years. 
%p - current proportion of population vaccinated.
%temp - tempory store of vaccinated and infectious individuals of a
%particular vaccination scheme.
function y = proportionidealBSVIRcost ()
costpermeaslescase = 307;
%This may be delivery cost as well as the MMR so 13.4+8.7=22.1
costpervaccintion = 22.1;
associatedcostpervaccination = 2.08;
prop = [0:1/100:1];
l = length(prop);
c = 1;
infections = zeros (1,l);
vaccinations = zeros (1,l);
%To save on time run the system for 10 years, then 20 years at 85%
%vaccination to determine the intial condiations once. Then these initial
%conditions can be used again in each model. (Avoiding repeat calculations)
intromassvaccyear = 10;
massvaccperiod = 20;
massvaccprop = 0.85;
br0 = 0;
bv0 = 0;
s0 = 6400000;
v0 = 0;
i0 = 600000;
r0 = 53000000;
u0 = [br0;bv0;s0;v0;i0;r0];
U = solveidealBSVIR (u0,massvaccprop, intromassvaccyear, (intromassvaccyear + massvaccperiod));
%We don't care about the cost we incur in the past so can ignore the past.
postmassvacc = U(:,end);
while (c<=l)
    p = prop(c);
    temp = solvenewinfectiousandvaccidealBSVIR(postmassvacc,p,0,50);
    infections(c)=sum(temp(1,:));
    vaccinations(c)=sum(temp(2,:));
    c = c+1;
end
totalinfectioncost = costpermeaslescase*infections;
totalvaccinationcost = (costpervaccintion+associatedcostpervaccination)*vaccinations;
totalcost = totalinfectioncost+totalvaccinationcost;
y = totalcost;
plot (prop,totalinfectioncost,'r',prop,totalvaccinationcost,'b',prop,totalcost,'k');
legend ('Infections','Vaccinations','Overall');
h_legend=legend ('Infections','Vaccinations','Overall');
set(h_legend,'FontSize',14);
x_label = xlabel('Proportion Vaccinated after mass-vaccination');
set(x_label,'FontSize',18);
y_label = ylabel('Total Cost US$ 2001');
set(y_label,'FontSize',16);