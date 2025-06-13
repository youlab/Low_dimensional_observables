%simulate target (focal) population dynamics using bounded gLV model
%three scenarios are simulated, within the same parameter set:

clear all; close all;
n_ori=100;%community size
model_index=1;
rng(model_index);
gamma=-rand(n_ori,n_ori)*0.4;
gamma=gamma - diag(diag(gamma));
alpha=rand(n_ori,n_ori)*0.4;
alpha=alpha - diag(diag(alpha));
sigma=rand(n_ori,1)*0.1+0.05;
mu=rand(n_ori,1)*0.5+0.5;
tspan=[0 300-6];%simulated time span to get steady state of the background population
y0=rand(n_ori,1)*0.2;
[t,results]=ode15s(@(t,y)bounded_gLV(t,y,gamma,alpha,sigma,mu),tspan,y0);

%%determine if the simulaiton is reasonable
if max(results,[],"all")>1
    display("simulation out of bound")
end
figure(1);
for i = 1:n_ori
    plot(t,results(:,i));hold on;
end
xlabel("time");

figure(2);
endpoint = results(end,:);
bar(endpoint);
survival = sum(endpoint>0.01);
fprintf("survival fraction: %.2f\n",survival/n_ori)

[~, sorted_idx] = sort(endpoint, 'descend');
gamma = gamma(sorted_idx,sorted_idx);
alpha = alpha(sorted_idx,sorted_idx);
sigma = sigma(sorted_idx);
mu = mu(sorted_idx);
y0 = y0(sorted_idx);

% check if the result is sorted correctly
[t,results]=ode15s(@(t,y)bounded_gLV(t,y,gamma,alpha,sigma,mu),tspan,y0);

figure(3);
for i = 1:n_ori
    plot(t,results(:,i));hold on;
end
xlabel("time");

figure(4);
endpoint = results(end,:);
bar(endpoint); hold on;
potential_target = sum(endpoint>0.1);

target_idx = randsample(1:potential_target, 10);
bar(target_idx,endpoint(target_idx));

figure(5);
for i = target_idx
    plot(t,results(:,i));hold on;
end
xlabel("time");

all_idx = 1:100;
background_idx = setdiff(all_idx, target_idx, 'stable');
% Concatenate: target first, then rest
final_order = [target_idx, background_idx];

gamma = gamma(final_order,final_order);
alpha = alpha(final_order,final_order);
sigma = sigma(final_order);
mu = mu(final_order);
y0 = y0(final_order);

% check if the result is sorted correctly
[t,results]=ode15s(@(t,y)bounded_gLV(t,y,gamma,alpha,sigma,mu),tspan,y0);

figure(6);
endpoint = results(end,:);
bar(endpoint); hold on;

figure(7);
for i = 1:10
    plot(t,results(:,i));hold on;
end
xlabel("time");

foldername = sprintf('./I%i', model_index);
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

filename = sprintf("./I%i/bounded_O%i_I%i.mat",model_index,n_ori,model_index);
% save(filename,"n_ori","gamma","alpha","sigma","mu","y0","-v7");
fprintf("Saved file: %s\n", fullfile(pwd, filename));