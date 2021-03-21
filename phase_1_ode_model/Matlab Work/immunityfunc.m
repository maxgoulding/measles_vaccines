%immunityfunc plots the immunity functions, i.e. the proportion of infants
%immune to measles as a function of time. You can vary the lambda the scale and n the shape
%values, for the gamma distribution of maternal antibodies.
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
%***x - the value required by the proportion function given interms of
%maternal antibodies outline in part 3***
function immunityfunc (n,lambda)
k = 1.1695;
r = 0.9216;
mu = 5.16;
T = 300;
%sato works around T = 200, similar to other PRN based tests
%lennon give HI titre threshold of 5
a = (0:1/48:1.25);
x = ((T.*exp(mu.*a)/(2^k)).^(1/r));
prop = 1- gamcdf(x,n,lambda);
%Plot it as months rather than years, as early vaccination scheduals are
%month based.
plot ((a*12),prop,'r');
axis ([0 15 0 1]);
axis square;
set(gca,'XTick',0:3:15);
set(gca,'YTick',0:0.2:1);
x_label = xlabel('Age (Months)');
set(x_label,'FontSize',14);
y_label = ylabel('Proportion immune');
set(y_label,'FontSize',14);