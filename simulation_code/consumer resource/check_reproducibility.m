% check_reproducibility.m
% Script to verify reproducibility of MCRM simulations

clear; clc;

%% Choose experiment index and simulation to test
index = 1;
test_sim = 5;   % which simulation to reproduce

%% Load saved parameters
load(sprintf('./I%i/params.mat',index),'params');

numSpecies = params.num_species;
numResources = params.num_resources;

%% Load saved initial conditions
load(sprintf('./I%i/inits.mat',index),'init');

%% Load saved outputs
species_file = sprintf('./I%i/mcrm_species.txt',index);
resource_file = sprintf('./I%i/mcrm_resource.txt',index);

species_saved = readmatrix(species_file);
resources_saved = readmatrix(resource_file);

%% Re-run selected simulation

params.x0 = init(test_sim,:)';     % restore initial condition
r = run_mcrm(params);              % run solver again

species_new = r.species;
resources_new = r.resources;

%% Extract original saved result

species_original = species_saved((test_sim-1)*numSpecies+1:test_sim*numSpecies,:)';
resources_original = resources_saved((test_sim-1)*numResources+1:test_sim*numResources,:)';

%% Compare results

species_error = max(abs(species_original - species_new), [], 'all');
resource_error = max(abs(resources_original - resources_new), [], 'all');

fprintf('\nReproducibility Check for simulation %d\n', test_sim);
fprintf('Max species difference   : %.6e\n', species_error);
fprintf('Max resource difference  : %.6e\n\n', resource_error);

%% Optional: visualize comparison

figure;

subplot(2,1,1)
plot(species_original,'b')
hold on
plot(species_new,'r--')
title('Species abundance comparison')
legend('Saved','Re-run')

subplot(2,1,2)
plot(resources_original,'b')
hold on
plot(resources_new,'r--')
title('Resource abundance comparison')
legend('Saved','Re-run')
