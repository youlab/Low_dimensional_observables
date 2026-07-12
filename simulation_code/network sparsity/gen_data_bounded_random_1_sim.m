close all; clear all;

% time step
% model 2: dt = 0.4
% model 3: dt = 0.4
% model 4: dt = 0.5
% model 5: dt = 0.5
% model 6: dt = 0.8
% model 7: dt = 2
% model 8: dt = 6

filename = "random_init.txt";
init_random = readmatrix(filename);
tspan=[0 300];
n_ori=100;
n_sim=10000;

index=4; dt=0.5;

fprintf("n_ori=%i, index=%i\n", n_ori, index);
load(sprintf("./I%i/bounded_O%i_I%i.mat",index,n_ori,index));

target10 = cell(n_sim, 1);
data_all = cell(n_sim, 1);
parfor i = 1:n_sim
    init=init_random(i,:)';
    sol=ode15s(@(t,y)bounded_gLV(t,y,gamma,alpha,sigma,mu),tspan,init);
    time=0:1:49;
    y=deval(sol,time*dt);
    target10{i} = y(1:10,:);
    data_all{i} = y;
end

target10 = vertcat(target10{:});
filename10=sprintf("./I%i/bgLV_T10_random.txt",index);
%writematrix(target10,filename10,'Delimiter','space');

data_all = vertcat(data_all{:});
filename=sprintf("./I%i/bgLV_random_all.txt",index);
writematrix(data_all,filename,'Delimiter','space');

figure(1);
x=0:1:49;
for i = 1:4
    subplot(2, 2, i);  % 2 rows, 2 columns, subplot i
    y = target10(10*(i-1)+1:10*i,:);
    for j = 1:10
        plot(x, y(j,:));hold on;
    end
end
