%Fits a binomial distribution linear regression to the data given. The
%data should contain the dates in which the sample was taken, the number of
%samples that we successfully immune (i.e. successes) and the number of
%samples that were taken.
%Date - The date in months when the samples were taken
%immune - the number of samples that we successfully tested as immune.
%sample - the number of samples taken.
%x(1) - Date in months
%x(2) - Rounded to 3 s.f. Linear Regression Binomial fit. Rounded to 3 s.f.
%x(3) - Lower bound of 95% Confidence interval. Rounded to 3 s.f. 
%x(4) - Upper bound of 95% Confidence interval. Rounded to 3 s.f. 
function x = binomialfit(date,sample,immune)
months = [0:1:15]';
[b,dev,stats] = glmfit (date, [immune sample], 'binomial', 'link', 'probit');
[immunefit,dlow,dhigh] = glmval(b, months,'probit',stats);
immunefitlow = immunefit-dlow;
immunefithigh = immunefit+dhigh;
plot(months,immunefit,'-',months,immunefitlow,'k--',months,immunefithigh,'k--');
axis ([0 15 0 1]);
axis square;
set(gca,'XTick',0:1:15);
set(gca,'YTick',0:0.2:1);
x_label = xlabel('Age (Months)');
set(x_label,'FontSize',18);
y_label = ylabel('Proportion immune');
set(y_label,'FontSize',18);
h_legend=legend('Linear Fit','95% Condidence');
set(h_legend,'FontSize',18);
roundedimmunefit = round(immunefit*1000)/1000;
roundedimmunefitlow = round(immunefitlow*1000)/1000;
roundedimmunefithigh = round(immunefithigh*1000)/1000;
x = [months,roundedimmunefit,roundedimmunefitlow,roundedimmunefithigh];