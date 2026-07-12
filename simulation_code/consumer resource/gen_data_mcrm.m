% MiCRM library based on: https://github.com/jgoldford/mcrm (J. Goldford)
% See ./README.md for attribution.
close all; clear all;
index = 3;
n_sim = 10000;
rng(index+43);
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

% run ODE solver for ecosystem over n_sim different initial conditions
init = zeros(n_sim, numSpecies+numResource);
init(:,1:numSpecies) = rand(n_sim, numSpecies)*1e5;
init(:, numSpecies+1:numSpecies+numResource) = rand(n_sim,numResource)*5;

% save parameters and initial conditions for reproducibility
save(sprintf('./I%i/params.mat',index),'params');
save(sprintf('./I%i/inits.mat',index),'init');

community_data = cell(n_sim, 1);
resource_data = cell(n_sim,1);
tic

for i = 1:n_sim
    
    params.x0 = init(i,:)';
    r = run_mcrm(params);
    
    community_data{i} = r.species';
    resource_data{i} = r.resources';

    % Progress report every 10 iterations (change if needed)
    if mod(i,100) == 0 || i == n_sim
        elapsed = toc;
        rate = i / elapsed;
        remaining = (n_sim - i) / rate;
        
        fprintf('Progress: %d/%d (%.1f%%) | Elapsed: %.1fs | Remaining: %.1fs\n',...
            i, n_sim, 100*i/n_sim, elapsed, remaining);
    end
end

figure()
for i = 1:3
    subplot(3,2,2*(i-1)+1);
    plot(community_data{i}');
    subplot(3,2,2*i);
    plot(resource_data{i}');
end

community_data = vertcat(community_data{:});
resource_data = vertcat(resource_data{:});


filename=sprintf("./I%i/mcrm_species.txt",index);
writematrix(community_data,filename,'Delimiter','space');

filename2=sprintf("./I%i/mcrm_resource.txt",index);
writematrix(resource_data,filename2,'Delimiter','space');

