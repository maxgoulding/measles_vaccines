function x = freqtopoint(interval,freq)
freqcopy = freq;
numinterval = length (interval);
i=1;
j=1;
while (i<=numinterval)
    if (freqcopy(i)~=0)
       x(j) = interval(i);
       freqcopy(i) = freqcopy(i)-1;
       j=j+1;
    else
       i=i+1; 
    end
end
