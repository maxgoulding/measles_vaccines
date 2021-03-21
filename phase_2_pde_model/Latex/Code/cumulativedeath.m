%Takes the average life expectancy within a population and take the age at
%which death starts occuring and returns the cumulative deaths within the
%population. This is a mixed distribution with uniform death before
%startage and exponential distribution afterwarss. Produces a graph.
%Returns values up to twice the avglife.
%deatth starts. Avglife must be greater than startage.
%startage - The age at which death starts occuring within the population.
%avglife - The average life expectancy within the population.
function [x,y] = cumulativedeath(startage,avglife)
remaining = avglife - startage;
x = [0:1:(avglife*2)];
xprime = [0:1:((2*avglife)-startage)];
yprime = expcdf(xprime,remaining);
yzeros = zeros(1,startage);
y = [yzeros,yprime];
plot(x,y,'-')
axis ([0 (2*avglife) 0 1]);
set(gca,'XTick',0:10:(2*avglife));
set(gca,'YTick',0:0.1:1);
x_label = xlabel('Age (Years)');
set(x_label,'FontSize',18);
y_label = ylabel('Proportion of cohort dead');
set(y_label,'FontSize',18);