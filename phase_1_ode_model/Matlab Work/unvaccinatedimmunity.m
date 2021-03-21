%Returns the proportion of infants with maternal immunity born from
%unvaccinated mothers after a certain period of time after birth in years.
function y = unvaccinatedimmunity(years)
n = 1.26;
lambda = 3075;
y = immunityyears(years,n,lambda);