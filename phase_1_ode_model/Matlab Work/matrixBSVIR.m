%Differential Model for realistic BVSIR system dynamics storing the
%differences as a matrix.
%u - matrix containing values
%totalI-Total in Infectious Group
%Population Assumptions
%N - total population size.
%lifexpt - Avg life expectancy. Years
%alpha - birth/death rate. Per person per year.
%infectiousperiod - Avg infectious period of measles.
%gamma - recovery rate. Per person per year.
%infectionage - avergage age of infection (without vaccination). years
%delta - Force of infection. per person person year
%beta - transmission rate. per person per year
%k - the proportion of people who make up the birth-able groups
%l - length of u
%c - counter
%Br-Birth Immunity Group. Of recovered mothers
%Bv-Birth Immunity Group. Of vaccinated mothers.
%S-Susceptable Group
%V-Vaccinated Group
%I-Infectious Group
%R-Recovered Group
function uprime = matrixBSVIR(t,u)
totalI=sum(u(5,:));
N = 60000000;
lifeexpt = 70;
alpha = 1/lifeexpt;
infectiousperiod = 12/365;
gamma = 1/infectiousperiod;
infectionage = 5;
delta = 1/infectionage;
beta = ((delta+alpha)*(alpha+gamma))/(alpha*N);
l = length(u);
c = 1;
brprime = zeros(1,l);
bvprime = zeros(1,l);
sprime = zeros(1,l);
vprime = zeros(1,l);
iprime = zeros(1,l);
rprime = zeros(1,l);
while (c<=l)
    Br = u(1,c);
    Bv = u(2,c);
    S = u(3,c);
    V = u(4,c);
    I = u(5,c);
    R = u(6,c);
    brprime(c) = -alpha*Br;
    bvprime(c) = -alpha*Bv;
    sprime(c) = -(beta*S*totalI)-(alpha*S);
    vprime(c) = -alpha*V;
    iprime(c) = (beta*S*totalI)-(gamma*I)-(alpha*I);
    rprime(c) = (gamma*I)-(alpha*R);
    c=c+1;
end
uprime = [brprime;bvprime;sprime;vprime;iprime;rprime];