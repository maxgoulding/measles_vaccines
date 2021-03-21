%Takes a set of maternal antibody data and plots the resulting immunity
%function as well as the lower and upper confidence bounds that would 
%be obtained from parameterising the data and inputting parameters into the
%immuniy function.
%T - Threshold of immunity
%data - Vector of datapoints of maternal antibody titre.
%est - estimated scale and shape parameters for the data
%conf - 95% lower and upper bounds for the scale and shape parameters for
%the data
%k - constant term to translate maternal into infant antibody titre.
%r - linear coefficient to translate maternal into infant antibody titre.
%mu - decay rate of infant antibody titre
%x - input values for gammacdf.
function datatoimmunityfunc (T,data)
%Parameterise the data
[est,conf] = gamfit(data);
k = 1.1695;
r = 0.9216;
mu = 0.43;
a = (0:0.25:15);
x = ((T.*exp(mu.*a)/(2^k)).^(1/r));
prop = 1- gamcdf(x,est(1),est(2));
proplow = 1- gamcdf(x,conf(1,1),conf(1,2));
prophigh = 1- gamcdf(x,conf(2,1),conf(2,2));
plot (a,prop,'-',a,proplow,':',a,prophigh,':');
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