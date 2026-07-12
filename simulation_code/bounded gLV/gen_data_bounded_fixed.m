close all; clear all;
tspan=[0 300];
n_ori=100;
n_sim=10000;
index=3;
fprintf("n_ori=%i, index=%i\n", n_ori, index);
load(sprintf("./I%i/bounded_O%i_I%i.mat",index,n_ori,index));
for n_target = [4,6,7,9]
    rng(10*n_target+index);
    n_background=n_ori-n_target;
    %%% fixed background simulation
    init_fixed=rand(n_sim,n_target)*0.2;
    target_all = cell(n_sim, 1);
    parfor i = 1:n_sim
        init = [init_fixed(i,:)'; y0(n_target+1:end)];
        sol=ode15s(@(t,y)bounded_gLV(t,y,gamma,alpha,sigma,mu),tspan,init);
        time=0:6:300-6;
        y=deval(sol,time);
        target_all{i} = y(1:n_target,:);
    end
    target_fixed = vertcat(target_all{:});
    size(target_fixed)
    filename=sprintf("./I%i/bgLV_B%i_T%i_fixed.txt",index,n_background,n_target);
    writematrix(target_fixed,filename,'Delimiter','space');
end