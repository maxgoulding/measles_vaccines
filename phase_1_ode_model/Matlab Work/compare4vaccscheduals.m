function y = compare4vaccscheduals()
currentvacc = proportionrealisticdoubleBSVIRcost(1,4);
ninemonthdouble = proportionrealisticdoubleBSVIRcost(273/365,4);
singlevacc12month = proportionrealisticBSVIRcost(1);
singlevacc9month = proportionrealisticBSVIRcost(273/365);
y = [currentvacc;ninemonthdouble;singlevacc12month;singlevacc9month];
t = [0:1/100:1];
plot (t,y(1,:),'r',t,y(2,:),'b',t,y(3,:),'k',t,y(4,:),'g');
h_legend=legend('Vaccination at 12 months and 4 years (Current)','Vaccination at 9 months and 4 years','Vaccination at 12 months','Vaccination at 9 months');
set(h_legend,'FontSize',14);
x_label = xlabel('Proportion Vaccinated after mass-vaccination');
set(x_label,'FontSize',18);
y_label = ylabel('Total Cost US$ 2001');
set(y_label,'FontSize',16);
%t = [0.8:1/100:1];
%plot (t,A(1,81:end),'r',t,A(2,81:end),'b',t,A(3,81:end),'k',t,A(4,81:end),'g');
%axis([0.8 1 0 4500000000])
%h_legend=legend('Vaccination at 12 months and 4 years (Current)','Vaccination at 9 months and 4 years','Vaccination at 12months','Vaccination at 9months');
%set(h_legend,'FontSize',14);
%x_label = xlabel('Proportion Vaccinated after mass-vaccination');
%set(x_label,'FontSize',14);
%y_label = ylabel('Total Cost US$ 2001');
%set(y_label,'FontSize',14);