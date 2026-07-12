close all; clear all;
tspan=[0 300];
n_ori=100;
n_sim=100000;
index=3;
fprintf("n_ori=%i, index=%i\n", n_ori, index);
load(sprintf("./I%i/dispersal_O%i_I%i.mat",index,n_ori,index));
rng(index);

%%% random background simulation
init_random=rand(n_sim,n_ori)*0.2;
filename=sprintf("./I%i/dgLV_random_init.txt",index);
writematrix(init_random,filename,'Delimiter','space');
target2 = cell(n_sim, 1);
target3 = cell(n_sim, 1);
target5 = cell(n_sim, 1);
target8 = cell(n_sim, 1);
parfor i = 1:n_sim
    init=init_random(i,:)';
    sol=ode15s(@(t,y)dispersal_gLV(t,y,gamma,mu,D),tspan,init);
    time=0:6:294;
    y=deval(sol,time);
    target2{i} = y(1:2,:);
    target3{i} = y(1:3,:);
    target5{i} = y(1:5,:);
    target8{i} = y(1:8,:);
end

target2 = vertcat(target2{:});
filename2=sprintf("./I%i/dgLV_B%i_T%i_random.txt",index,98,2);
writematrix(target2,filename2,'Delimiter','space');

target3 = vertcat(target3{:});
filename3=sprintf("./I%i/dgLV_B%i_T%i_random.txt",index,97,3);
writematrix(target3,filename3,'Delimiter','space');

target5 = vertcat(target5{:});
filename5=sprintf("./I%i/dgLV_B%i_T%i_random.txt",index,95,5);
writematrix(target5,filename5,'Delimiter','space');

target8 = vertcat(target8{:});
filename8=sprintf("./I%i/dgLV_B%i_T%i_random.txt",index,92,8);
writematrix(target8,filename8,'Delimiter','space');
