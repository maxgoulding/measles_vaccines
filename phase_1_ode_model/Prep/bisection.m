function [root iteration] = bisection (f,a,b,acc)
    fin = inline(char(f));
    c = (a+b)/2;
    count = 0;
    fc = feval (fin,c);
    while ((fc>acc)||(fc<-acc))
       if (fc>0)
           %b=c in this case
           b=c;
       else
           %a=c in this case
           a=c;
       end
       c=(a+b)/2;
       fc = feval (fin,c);
       count = count + 1;
    end
    root = c;
    iteration = count;