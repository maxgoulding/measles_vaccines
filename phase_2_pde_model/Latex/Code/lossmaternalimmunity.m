%The function takes the proporion of those indiviudals with maternal
%immunity after a period of time and translates this into the proportion
%that maintain maternal immunity at each stage. 
%date - Vector of dates/ages
%prop - The proportion of individuals with maternal immunity at date.
function x =lossmaternalimmunity (date,prop)
l=length(date);
y = zeros(l,1);
y(1)=prop(1);
i=2;
while (i<l)&&(prop(i-1)~=0)
    y(i) = prop(i)/prop(i-1);
    i=i+1;
end
x=[date,y];