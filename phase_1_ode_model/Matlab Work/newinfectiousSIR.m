%Differential Model estimating the new infectious individuals in the simple
%SIR model.
%S-Susceptable Group
%I-Infectious Group
%R-Recovered Group
%The Population assumptions
%lifexpt - Avg life expectancy. Years
%infectiousperiod - Avg infectious period of measles. Years
%infectionage - avergage age of infection (without vaccination). Years
%alpha - birth/death rate. Per person per year.
%gamma - recovery rate. Per person per year.
%delta - Force of infection. per person person year
%beta - transmission rate. per person per year
function uprime = newinfectiousSIR(t,u)
S=u(1);
I=u(2);
R=u(3);
N = S+I+R;
lifeexpt = 70;
alpha = 1/lifeexpt;
infectiousperiod = 12/365;
gamma = 1/infectiousperiod;
infectionage = 5;
delta = 1/infectionage;
beta = ((delta+alpha)*(alpha+gamma))/(alpha*N);
sprime = 0;
iprime = (beta*S*I);
rprime = 0;
uprime = [sprime;iprime;rprime];