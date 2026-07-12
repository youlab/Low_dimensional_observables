% MiCRM library based on: https://github.com/jgoldford/mcrm (J. Goldford)
% See ./README.md for attribution.
close all; clear all;
index = 1;
rng(index+45);
% define hyperparamters
% the total number of species in your "world"

Total_Number_Of_Species = 1000;

% the number of species in your ecosystem (i.e. inocula)
numSpecies = 100;

% the total number of metabolic families
numGroups = 10;

% the total number of resources (must equal or exceed numGroups)
numResource = 30;

% dirichlet hypoerparameter (controls the variablity of the family
% metabolism)
Total = 100;

% get a secretion matrix for everyone
D = getMetabolism(numResource,'rand',1/numResource);

% determine how much of the community you want to be specialists vs.
% generalists
% more generalist, set to 1/M, and more specialist set to 0.9
specialist = 0.7;

specialistVariation = 0.2;

priors = arrayfun(@(x) getConsumerPriors(20,x,specialist,specialistVariation,numResource)',1:numGroups,'uni',0);
[out,~]=makePhyloConsumers(Total_Number_Of_Species,numResource,numGroups,priors,Total,true);

% if you don't want community struture, use the function
% getRandomConsumerMatrix


% randomly sample a sub population:
k = randsample(1:Total_Number_Of_Species,numSpecies);
% construct an ecosystem
params = mcrm_params(1,numSpecies,out.C(k,:),D,'eye','');

% run ODE solver for ecosystem
r = run_mcrm(params);

% endpoint abundances
endpoint_abund = r.species(end, :);

% sort species by endpoint abundance (descending)
[sorted_abund, sort_idx] = sort(endpoint_abund, 'descend');
species_timeseries = r.species(:, sort_idx);

% plot original time series
figure();
subplot(1,2,1);
plot(species_timeseries(:,1:10));
ylabel('species abundance');
xlabel('time');

subplot(1,2,2);
plot(r.resources);
ylabel('resource abundance');
xlabel('time');

% plot sorted endpoint abundances
figure();
subplot(1,2,1);
bar(sorted_abund(1:20));
xlabel('species rank (sorted by endpoint abundance)');
ylabel('endpoint abundance');
title('Sorted endpoint abundances');

% plot distribution of endpoint abundances
subplot(1,2,2);
histogram(endpoint_abund, 30);
xlabel('endpoint abundance');
ylabel('number of species');
title('Distribution of endpoint abundances');
