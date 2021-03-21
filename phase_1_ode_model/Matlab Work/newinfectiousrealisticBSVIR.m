%Differential Model estimating the new infectious individuals in the
%realistic BSVIR model.
%Br-Birth Immunity Group. Of recovered mothers
%Bv-Birth Immunity Group. Of vaccinated mothers.
%S-Susceptable Group
%V-Vaccinated Group
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
%k - the proportion of people who make up the birth-able groups
function uprime = newinfectiousrealisticBSVIR(t,u)
Br=u(1);
Bv=u(2);
S=u(3);
V=u(4);
I=u(5);
R=u(6);
N = Br+Bv+S+V+I+R;
lifeexpt = 70;
alpha = 1/lifeexpt;
infectiousperiod = 12/365;
gamma = 1/infectiousperiod;
infectionage = 5;
delta = 1/infectionage;
beta = ((delta+alpha)*(alpha+gamma))/(alpha*N);
brprime = 0;
bvprime = 0;
sprime = 0;
vprime = 0;
iprime = (beta*S*I);
rprime = 0;
uprime = [brprime;bvprime;sprime;vprime;iprime;rprime];