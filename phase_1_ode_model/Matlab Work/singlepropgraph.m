% Uses glmfit to fit a log linear model to the data provided and draw the
% appropriate graph for it. It also labels the graph appropriately. It is
% equalant to doublepropgraph except it works with only one set of data.
% data - is matrix with width 3. Column 1 represents x, column 2 y,
% column 3 is n.
% x - is the intervals at which the readings were taken
% yu - is the number of successes from infants of unvaccinated mothers
% nu - is the number of trials from of unvaccinated mothers
function singlepropgraph (data)
x = data(:,1);
n = data(:,2);
y = data(:,3);
x = x(~isnan(n));
y = y(~isnan(n));
n = n(~isnan(n));
b= glmfit(x, [y,n], 'binomial', 'link', 'logit');
yfit = glmval(b,x,'logit','size', n);
plot (x,yfit./n,'r-',x,y./n,'r.');
axis ([0 15 0 1]);
axis square;
set(gca,'XTick',0:3:15);
legend ('Infants of Unvaccinated Mothers');
xlabel('Age (Months)')
ylabel('Proportion immune')