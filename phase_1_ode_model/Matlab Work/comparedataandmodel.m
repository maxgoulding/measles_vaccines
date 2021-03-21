%Plots a graph that contains a scatter plot obtained from experimental data
%with a curve given with the modelled immunity function given the
%parameters n and lambda.

%Parameters
%data - n x 5 matrix. 1st column is age (months).2nd column is number of
%infants of unvaccinated mothers. 3rd column is number of infants of
%unvaccinated mothers with immunity. 4th column is number of
%infants of vaccinated mothers. 5th column is number of infants of
%vaccinated mothers with immunity.
%%T - Threshold for immunity
%lambda1 - Scale parameter of first gamma distribution. (Red)
%n1 - Shape parameter of first gamma distribution. (Red)
%lambda2 - Scale parameter of second gamma distribution. (Blue)
%n2 - Shape parameter of second gamma distribution. (Blue)

%Variables used for scatter plot
%ae - age vector. The ages at which data was taken.
%nu - number of infants of unvaccinated mothers vector.
%iu - number of infants of unvaccinated mothers with immunity vector.
%nv - number of infants of vaccinated mothers vector.
%iv - number of infants of vaccinated mothers with immunity vector.

%Variables used for modelled immunity function
%r - ratio of linear relationship between log(cord) and log(maternal).
%Taken as the corrected values from Transplacental transfer of measles and total IgG.
%k - constant term of the linear relationship between log(cord) and
%log(maternal). %Taken as the corrected values from Transplacental transfer
%of measles and total IgG.
%mu - rate of decay
%a - the time/age in months
%x - the value required by the proportion function given interms of
%maternal antibodies

function comparedataandmodel(data,T,n1,lambda1,n2,lambda2)
%Data for scatter plot
ae = data (:,1);
nu = data (:,2);
iu = data (:,3);
nv = data (:,4);
iv = data (:,5);
%Modelled Immunity Function
k = 1.1695;
r = 0.9216;
mu = 0.35;
a = (0:0.25:15);
x = ((T.*exp(mu.*a)/(2^k)).^(1/r));
prop1 = 1- gamcdf(x,n1,lambda1);
prop2 = 1- gamcdf(x,n2,lambda2);
%Plotting the graph
plot (a,prop1,'r',a,prop2,'b',ae,iu./nu,'rx',ae,iv./nv,'bx');
hlegend = legend ('Model 1', 'Model 2', 'Infants of Unvaccinated Mothers','Infants of Vaccinated Mothers');
set(hlegend,'FontSize',14);
axis ([0 15 0 1]);
axis square;
set(gca,'XTick',0:3:15);
set(gca,'YTick',0:0.2:1);
xlabel('Age (Months)')
ylabel('Proportion immune')
x_label = xlabel('Age (Months)');
set(x_label,'FontSize',14);
y_label = ylabel('Proportion immune');
set(y_label,'FontSize',14);