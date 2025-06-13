function y=dispersal_gLV(t,x,gamma,mu,D)
x(x<0)=0;
negative_interactions=gamma*x;
y=mu.*x.*(1-x-negative_interactions)+D;
end