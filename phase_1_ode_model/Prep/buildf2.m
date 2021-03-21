function f2 = buildf2(N)
x =0:1/N:1;
y=x;
f = max(x,1-x)'*max(y,1-y);
f2 = reshape(f(2:N,2:N),(N-1)^2,1);