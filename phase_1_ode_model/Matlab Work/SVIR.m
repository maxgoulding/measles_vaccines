%Differential Model for SVIR
%S-Susceptable Group
%V-Vaccinated Group
%I-Infectious Group
%R-Recovered Group
%lifexpt - Avg life expectancy. Years
%alpha - birth/death rate. Per person per year.
%infectiousperiod - Avg infectious period of measles.
%gamma - recovery rate. Per person per year.
%infectionage - avergage age of infection (without vaccination). years
%delta - Force of infection. per person person year
%beta - transmission rate. per person per year
%p - proportion of newborns vaccinated
function uprime = SVIR(p,t,u)
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
sprime = ((1-p)*alpha*(N))-(alpha*S)-(beta*S*I);
vprime = (p*alpha*(N))-(alpha*V);
iprime = (beta*S*I)-(gamma*I)-(alpha*I);
rprime = (gamma*I)-(alpha*R);
uprime = [sprime;vprime;iprime;rprime];