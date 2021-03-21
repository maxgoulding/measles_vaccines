%Returns the proportion of infants with maternal immunity born from
%vaccinated mothers after a certain period of time after birth in years.
function y = vaccinatedimmunity(years)
n = 1.26;
%This is a guess and must be justified.
lambda = 1000;
y = immunityyears(years,n,lambda);