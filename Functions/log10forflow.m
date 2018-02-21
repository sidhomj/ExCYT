function y=log10forflow(x)

xneg=x<=0;
xmod=subplus(x);
xmod=log10(xmod);
bottom=min(xmod(xmod~=-Inf));
xmod(xneg)=bottom;
y=xmod;

end