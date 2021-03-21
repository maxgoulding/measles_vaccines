function f1 = buildf1(N)
x =0:1/N:1;
y=x;
f = sin(pi*x)'*sin(pi*y);
f1 = reshape(f(2:N,2:N),(N-1)^2,1);