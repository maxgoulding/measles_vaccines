%The function takes data containing the age of sampling, number of infants 
%and number with immunity from both vaccinated and unvaccinated 
%mothers, in the form of a matrix, and returns scatter graph of the
%proportion of infants immune as a function of their age.
%data - n x 5 matrix. 1st column is age (months).2nd column is number of
%infants of unvaccinated mothers. 3rd column is number of infants of
%unvaccinated mothers with immunity. 4th column is number of
%infants of vaccinated mothers. 5th column is number of infants of
%vaccinated mothers with immunity.
%a - age vector. The ages at which data was taken.
%nu - number of infants of unvaccinated mothers vector.
%iu - number of infants of unvaccinated mothers with immunity vector.
%nv - number of infants of vaccinated mothers vector.
%iv - number of infants of vaccinated mothers with immunity vector.
function datatoscatter (data)
a = data (:,1);
nu = data (:,2);
iu = data (:,3);
nv = data (:,4);
iv = data (:,5);
plot (a,iu./nu,'rx',a,iv./nv,'bx');
axis ([0 15 0 1]);
axis square;
set(gca,'XTick',0:3:15);
set(gca,'YTick',0:0.2:1);
hlegend = legend ('Infants of Unvaccinated Mothers','Infants of Vaccinated Mothers');
set(hlegend,'FontSize',14);
x_label = xlabel('Age (Months)');
set(x_label,'FontSize',14);
y_label = ylabel('Proportion immune');
set(y_label,'FontSize',14);