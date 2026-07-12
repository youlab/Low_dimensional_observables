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
model_index=4;
load(sprintf("./I%i/bounded_O%i_I%i.mat",model_index,n_ori,model_index));

dt_values = [0.4,0.4,0.5,0.5,0.8,2,6];
dt = dt_values(model_index-1);
p_values = [0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5];
p = p_values(model_index-1);
% Convert adjacency matrix to graph object
G = graph(network);
figure(1);
plot(G);
title(sprintf('Random network, p = %.3f', p));

% Degree of each node
node_id = (1:n_ori)';

% If network is binary adjacency:
degree_values = degree(G);   % same as sum(network, 2) for undirected binary networks

% Save node ID and degree to txt file
degree_data = [node_id, degree_values];

filename = sprintf("./I%i/network_degrees_I%i.txt", model_index, model_index);
writematrix(degree_data, filename, 'Delimiter', 'tab');
figure(2);
histogram(degree_values, 'BinMethod', 'integers');
xlabel('Degree');
ylabel('Number of nodes');
title(sprintf('Degree distribution, p = %.3f', p));

% Closeness centrality of each node
closeness_values = centrality(G, 'closeness');
% Save node ID and centrality to txt file
closeness_data = [node_id, closeness_values];
filename = sprintf("./I%i/closeness_centrality_I%i.txt", model_index, model_index);
writematrix(closeness_data, filename, 'Delimiter', 'tab');

figure(3);
histogram(closeness_values);
xlabel('Closeness centrality');
ylabel('Number of nodes');
title(sprintf('Closeness centrality distribution, p = %.3f', p));


% Connected components
component_id = conncomp(G);
component_sizes = histcounts(component_id, 1:(max(component_id)+1));

n_components = max(component_id);
largest_component_size = max(component_sizes);
largest_component_fraction = largest_component_size / n_ori;

fprintf('Number of connected components: %d\n', n_components);
fprintf('Largest component size: %d nodes\n', largest_component_size);
fprintf('Largest component fraction: %.2f\n', largest_component_fraction);

figure(4);

histogram(component_sizes, 'BinMethod', 'integers');

xlabel('Component size');
ylabel('Number of components');
title(sprintf('Connected component sizes, p = %.3f', p));

gamma = gamma .* network;
alpha = alpha .* network;
tspan=[0 300-6];%simulated time span to get steady state of the background population
y0=rand(n_ori,1)*0.2;
sol=ode15s(@(t,y)bounded_gLV(t,y,gamma,alpha,sigma,mu),tspan,y0);
time=0:1:49; time = time*dt;
y=deval(sol,time);
%%determine if the simulaiton is reasonable
if max(y,[],"all")>1
    display("simulation out of bound")
end
figure(5);
for i = 1:n_ori
    plot(time,y(i,:));hold on;
end
xlabel("time");