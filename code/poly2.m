function [params,rsq,curve]=poly2(x,y)
fo=fitoptions('Method','NonLinearLeastSquares','lower',[-Inf,-Inf,-Inf],'upper',[Inf,Inf,Inf],'StartPoint',[1,1,1]);
ft=fittype('(a*(x^2))+(b*(x^1))+c','options',fo);
[curve,gof] = fit(x,y,ft);
a=curve.a;
b=curve.b;
c=curve.c;
params=[a;b;c];
rsq=gof.rsquare;
end