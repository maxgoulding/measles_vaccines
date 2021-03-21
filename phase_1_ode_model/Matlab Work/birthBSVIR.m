%Differential Model for realistic BVSIR births as a vector.
%Br-Birth Immunity Group. Of recovered mothers
%Bv-Birth Immunity Group. Of vaccinated mothers.
%S-Susceptable Group
%V- Vaccinated Group
%I-Infectious Group
%R-Recovered Group
%N - Total population size
%lifexpt - Avg life expectancy. Years
%alpha - birth/death rate. Per person per year.
%k - the proportion of people who make up the birth-able groups
function uprime = birthBSVIR(t,u)
Br=u(1);
Bv=u(2);
S=u(3);
V=u(4);
I=u(5);
R=u(6);
N = 60000000;
lifeexpt = 70;
alpha = 1/lifeexpt;
k = N/(R+S+V);
ui = unvaccinatedimmunity(0);
vi = vaccinatedimmunity(0);
brprime = ui*(alpha*k*R);
bvprime = vi*(alpha*k*V);
sprime = (alpha*k*S)+((1-vi)*(alpha*k*V))+((1-ui)*(alpha*k*R));
vprime = 0;
iprime = 0;
rprime = 0;
uprime = [brprime;bvprime;sprime;vprime;iprime;rprime];