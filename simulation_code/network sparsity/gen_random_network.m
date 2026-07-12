% Generate random interaction networks for bounded gLV model
% start off with I1 in the original fully connected network simulations
% model 2: p = 0;
% model 3: p = 0.01;
% model 4: p = 0.02;
% model 5: p = 0.05;
% model 6: p = 0.1;
% model 7: p = 0.2;
% model 8: p = 0.5;

clear all; close all;
rng(42);
n_ori=100; %community size
load(sprintf("./I1/bounded_O%i_I1.mat",n_ori));

model_index=8;
p_values = [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5];
p = p_values(model_index-1);
network = random_network(n_ori,p);
% Convert adjacency matrix to graph object
G = graph(network);
figure(1);
plot(G);
title(sprintf('Random network, p = %.3f', p));

% Degree of each node
degree = sum(network, 2);
figure(2);
histogram(degree, 'BinMethod', 'integers');
xlabel('Degree');
ylabel('Number of nodes');
title(sprintf('Degree distribution, p = %.3f', p));

% Connected components
component_id = conncomp(G);
component_sizes = histcounts(component_id, 1:(max(component_id)+1));

n_components = max(component_id);
largest_component_size = max(component_sizes);
largest_component_fraction = largest_component_size / n_ori;

fprintf('Number of connected components: %d\n', n_components);
fprintf('Largest component size: %d nodes\n', largest_component_size);
fprintf('Largest component fraction: %.2f\n', largest_component_fraction);

figure(3);

histogram(component_sizes, 'BinMethod', 'integers');

xlabel('Component size');
ylabel('Number of components');
title(sprintf('Connected component sizes, p = %.3f', p));

gamma = gamma .* network;
alpha = alpha .* network;
tspan=[0 300-6];%simulated time span to get steady state of the background population
y0=rand(n_ori,1)*0.2;
[t,results]=ode15s(@(t,y)bounded_gLV(t,y,gamma,alpha,sigma,mu),tspan,y0);

%%determine if the simulaiton is reasonable
if max(results,[],"all")>1
    display("simulation out of bound")
end
figure(4);
for i = 1:n_ori
    plot(t,results(:,i));hold on;
end
xlabel("time");

foldername = sprintf('./I%i', model_index);
if ~exist(foldername, 'dir')
    mkdir(foldername);
end

filename = sprintf("./I%i/bounded_O%i_I%i.mat",model_index,n_ori,model_index);
save(filename,"n_ori","network","gamma","alpha","sigma","mu","y0","-v7");
fprintf("Saved file: %s\n", fullfile(pwd, filename));

function A = random_network(n, p)
    A = triu(rand(n) < p, 1);
    A = A + A';
end