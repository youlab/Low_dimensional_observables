%simulate target (focal) population dynamics using bounded gLV model
%three scenarios are simulated, within the same parameter set:

clear all; close all;
n_ori=10;%community size
model_index=6;
rng(model_index);
gamma=-rand(n_ori,n_ori)*0.4;
gamma=gamma - diag(diag(gamma));
alpha=rand(n_ori,n_ori)*0.4;
alpha=alpha - diag(diag(alpha));
sigma=rand(n_ori,1)*0.1+0.05;
mu=rand(n_ori,1)*0.5+0.5;
tspan=[0 49];%simulated time span to get steady state of the background population
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

foldername = sprintf('./I%i', model_index);
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

filename = sprintf("./I%i/bounded_O%i_I%i.mat",model_index,n_ori,model_index);
save(filename,"n_ori","gamma","alpha","sigma","mu","y0","-v7");
fprintf("Saved file: %s\n", fullfile(pwd, filename));