close all; clear all;
tspan=[0 300];
n_ori=100;
n_sim=10000;
index=3;
fprintf("n_ori=%i, index=%i\n", n_ori, index);
load(sprintf("./I%i/bounded_O%i_I%i.mat",index,n_ori,index));
for n_target = [1,4,6,7,9]
    rng(10*n_target+index);
    n_background=n_ori-n_target;
    %%% random background simulation
    init_random=rand(n_sim,n_ori)*0.2;
    target_all = cell(n_sim, 1);
    parfor i = 1:n_sim
        init=init_random(i,:)';
        sol=ode15s(@(t,y)bounded_gLV(t,y,gamma,alpha,sigma,mu),tspan,init);
        time=0:6:300-6;
        y=deval(sol,time);
        target_all{i} = y(1:n_target,:);
    end
    target_random = vertcat(target_all{:});
    size(target_random)
    filename=sprintf("./I%i/bgLV_B%i_T%i_random.txt",index,n_background,n_target);
    writematrix(target_random,filename,'Delimiter','space');
    filename2=sprintf("./I%i/bgLV_B%i_T%i_random_init.txt",index,n_background,n_target);
    writematrix(init_random,filename2,'Delimiter','space');
end