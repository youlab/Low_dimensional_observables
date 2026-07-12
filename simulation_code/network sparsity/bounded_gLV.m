function y=bounded_gLV(t,x,gamma,alpha,sigma,mu)
x(x<0)=0;
positive_interactions=alpha*x;
negative_interactions=gamma*x;
death=sigma./(1+positive_interactions);
death(death<0)=0;
y=mu.*x.*(1-x-death+negative_interactions);
end