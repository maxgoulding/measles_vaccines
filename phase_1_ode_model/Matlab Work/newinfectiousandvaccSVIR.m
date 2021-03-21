%Differential Model estimating the new infectious and vaccinated 
%individuals in the SVIR model.
%S-Susceptable Group
%V-Vaccinated Group
%I-Infectious Group
%R-Recovered Group
%p - The proportion of newborns vaccinated
%The Population assumptions
%lifexpt - Avg life expectancy. Years
%infectiousperiod - Avg infectious period of measles. Years
%infectionage - avergage age of infection (without vaccination). Years
%alpha - birth/death rate. Per person per year.
%gamma - recovery rate. Per person per year.
%delta - Force of infection. per person person year
%beta - transmission rate. per person per year
function uprime = newinfectiousandvaccSVIR(p,t,u)
S=u(1);
V=u(2);
I=u(3);
R=u(4);
N = S+V+I+R;
lifeexpt = 70;
alpha = 1/lifeexpt;
infectiousperiod = 12/365;
gamma = 1/infectiousperiod;
infectionage = 5;
delta = 1/infectionage;
beta = ((delta+alpha)*(alpha+gamma))/(alpha*N);
sprime = 0;
vprime = (p*alpha*N);
iprime = (beta*S*I);
rprime = 0;
uprime = [sprime;vprime;iprime;rprime];