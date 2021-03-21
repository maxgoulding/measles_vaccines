%Differential Model for ideal BVSIR.
%p - proportion of newly susceptible infants vaccinated as they become
%susceptible.
%Br-Birth Immunity Group. Of recovered mothers
%Bv-Birth Immunity Group. Of vaccinated mothers.
%S-Susceptable Group
%V-Vaccinated Group
%I-Infectious Group
%R-Recovered Group
%Population Assumptions
%lossrecoveredmatimmunityage - Avg. loss of immunity of maternally derived
%antibodies for infants of recovered mother. Years.
%sigma - recovered exposure rate. Per person per year. 
%lossvaccinatedmatimmunityage - Avg. loss of immunity of maternally derived
%antibodies for infants of vaccinated mother. Years.
%xi - vaccinated exposure rate. Per person per year.
%lifexpt - Avg life expectancy. Years
%alpha - birth/death rate. Per person per year.
%infectiousperiod - Avg infectious period of measles.
%gamma - recovery rate. Per person per year.
%infectionage - avergage age of infection (without vaccination). years
%delta - Force of infection. per person person year
%beta - transmission rate. per person per year
%k - the proportion of people who make up the birth-able groups
function uprime = idealBSVIR(p,t,u)
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
brprime = (alpha*k*R)-(sigma*Br)-(alpha*Br);
bvprime = (alpha*k*V)-(xi*Bv)-(alpha*Bv);
sprime = ((1-p)*((sigma*Br)+(xi*Bv)))+((1-p)*alpha*k*S)-(beta*S*I)-(alpha*S);
vprime = (p*((sigma*Br)+(xi*Bv)+(alpha*k*S)))-(alpha*V);
iprime = (beta*S*I)-(gamma*I)-(alpha*I);
rprime = (gamma*I)-(alpha*R);
uprime = [brprime;bvprime;sprime;vprime;iprime;rprime];