%immunitydays return the proportion of infants still immune after a certain
%number of days given the day, scale and shape value for the gamma distn
%of the maternal antibody titre.
%Maternal antibodies are assumes to be linearly related by
%log(cord)=r*log(maternal)*k
%lambda - Scale parameter of gamma distribution.
%n - Shape parameter of gamma distribution.
%r - ratio of linear relationship between log(cord) and log(maternal).
%Taken as the corrected values from Transplacental transfer of measles and total IgG.
%k - constant term of the linear relationship between log(cord) and
%log(maternal). %Taken as the corrected values from Transplacental transfer
%of measles and total IgG.
%mu - rate of decay
%T - Threshold for immunity. Taken from standard mIU/ml
%a - the time/age in months
function y=immunityyears (years,n,lambda)
k = 1.1695;
r = 0.9216;
mu = 5.16;
T = 120;
x = ((T*exp(mu*years)/(2^k))^(1/r));
y = 1- gamcdf(x,n,lambda);
