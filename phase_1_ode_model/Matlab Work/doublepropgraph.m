% Uses glmfit to fit a log linear model to the data provided and draw the
% appropriate graph for it. Have set the default to title the graph and
% also to label the axis to proportion immune against age in months.
% data - is matrix with width 5. Column 1 represents x, column 2 yu,
% column 3 is nu, column 4 is yv and column 5 is nv.
% x - is the intervals at which the readings were taken
% yu - is the number of successes from infants of unvaccinated mothers
% nu - is the number of trials from of unvaccinated mothers
% yv - is the number of successes from infants of vaccinated mothers
% nv - is the number of trials from of vaccinated mothers
function doublepropgraph (data)
x = data(:,1);
nu = data(:,2);
yu = data(:,3);
nv = data(:,4);
yv = data(:,5);
xu = x(~isnan(nu));
xv = x(~isnan(nv));
yu = yu(~isnan(nu));
yv = yv(~isnan(nv));
nu = nu(~isnan(nu));
nv = nv(~isnan(nv));
bu = glmfit(xu, [yu,nu], 'binomial', 'link', 'logit');
bv = glmfit(xv, [yv,nv], 'binomial', 'link', 'logit');
yufit = glmval(bu,xu,'logit','size', nu);
yvfit = glmval(bv,xv,'logit','size', nv);
plot (xu,yufit./nu,'r-',xv,yvfit./nv,'b-',xu,yu./nu,'r.',xv,yv./nv,'b.');
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