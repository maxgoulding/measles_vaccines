%Differential Model estimating the new infectious and vaccinated 
%individuals in the ideal BSVIR model.
%Br-Birth Immunity Group. Of recovered mothers
%Bv-Birth Immunity Group. Of vaccinated mothers.
%S-Susceptable Group
%V-Vaccinated Group
%I-Infectious Group
%R-Recovered Group
%The Population assumptions
%lossrecoveredmatimmunityage - Avg. loss of immunity of maternally derived
%antibodies for infants of recovered mother. Years.
%sigma - recovered exposure rate. Per person per year. 
%lossvaccinatedmatimmunityage - Avg. loss of immunity of maternally derived
%antibodies for infants of vaccinated mother. Years.
%xi - vaccinated exposure rate. Per person per year.
%lifexpt - Avg life expectancy. Years
%infectiousperiod - Avg infectious period of measles. Years
%infectionage - avergage age of infection (without vaccination). Years
%alpha - birth/death rate. Per person per year.
%gamma - recovery rate. Per person per year.
%delta - Force of infection. per person person year
%beta - transmission rate. per person per year
%k - the proportion of people who make up the birth-able groups
%p - proportion of newly susceptible infants vaccinated as they become
%susceptible.
function uprime = newinfectiousandvaccidealBSVIR(p,t,u)
Br=u(1);
Bv=u(2);
S=u(3);
V=u(4);
I=u(5);
R=u(6);
N = Br+Bv+S+V+I+R;
lossrecoveredmatimmunityage = 1;
sigma = 1/lossrecoveredmatimmunityage;
lossvaccinatedmatimmunityage = 0.75;
xi = 1/lossvaccinatedmatimmunityage;
lifeexpt = 70;
alpha = 1/lifeexpt;
infectiousperiod = 12/365;
gamma = 1/infectiousperiod;
infectionage = 5;
delta = 1/infectionage;
beta = ((delta+alpha)*(alpha+gamma))/(alpha*N);
k = N/(R+S+V);
brprime = 0;
bvprime = 0;
sprime = 0;
vprime = (p*((sigma*Br)+(xi*Bv)+(alpha*k*S)));
iprime = (beta*S*I);
rprime = 0;
uprime = [brprime;bvprime;sprime;vprime;iprime;rprime];