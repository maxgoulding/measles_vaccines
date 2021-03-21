%Returns the proportion of infants who have lost maternal immunity between 
%two time intervals given they are born to unvaccinated mothers . Years.
%t1 - lower time interval
%t2 - upper time interval
%p1 - lower time interval immunity function
%p2 - upper time interval immunity function
function y = lossofunvaccinatedimmunity(t1,t2)
n = 1.26;
lambda = 3075;
p1 = immunityyears(t1,n,lambda);
p2 = immunityyears(t2,n,lambda);
%Numerical errors occur if we try to divide by 0 or close to 0. So we want
%everyone to lose it
if (p1 == 0)
    y=1;
else
    y =1-(p2/p1);
end