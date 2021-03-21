%Plots a graph that contains a scatter plot obtained from experimental data
%with a curve given with the modelled immunity function given experimental
%values for maternal antibody titre.

%Parameters
%modeldata - Vector of datapoints of maternal antibody titre.
%scatterdata - n x 3 matrix. 1st column is age (months).2nd column is number of
%infants. 3rd column is number of infantswith immunity.
%T - Threshold for immunity

%Variables used for scatter plot
%ae - age vector. The ages at which data was taken.
%n - number of infants vector.
%i - number of infants with immunity vector.

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
function modelscattercompare(T,modeldata,scatterdata)
%Data for scatter plot
ae = scatterdata (:,1);
n = scatterdata (:,2);
i = scatterdata (:,3);
%Parameterise the data
[est,conf] = gamfit(modeldata);
%Modelled Immunity Function
k = 1.1695;
r = 0.9216;
mu = 0.43;
a = (0:0.25:15);
x = ((T.*exp(mu.*a)/(2^k)).^(1/r));
prop = 1- gamcdf(x,est(1),est(2));
proplow = 1- gamcdf(x,conf(1,1),conf(1,2));
prophigh = 1- gamcdf(x,conf(2,1),conf(2,2));
%Plotting the graph
plot (a,prop,'-',a,proplow,':',a,prophigh,':',ae,i./n,'x');
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