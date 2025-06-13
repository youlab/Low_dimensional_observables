% Author: Emrah Simsek
% Latest update: 16-Feb-2024

% this code is to simulate growth dynamics on a n1-by-n2 cluster community.
% NNat # Background species randomly dispersed on the grid
% Nf # Target species also dispersed on a certain fraction of the grid sites
% This is a spatia:l gLV ODE model, a modified version
% of the model used in Feilun' s recent paper by an explicit consideretion
% of spatial configuration as a scaling factor of the interactions.
% Interactions within the same site is scaled by "1" (i.e. unscaled). So,
% their effective scaling could be the same as the nearest neighbor scaling
% if the nearest neighbor distance (dmin) is also "1".

% ix = (col - 1) * n + row; %linear index

% this code builds upon my "sgLV_v03b2.m" implementation by using a
% built-in ode solver, instead of using manual Euler method.
% here I use Zhengqing' s code snippet to write data every certain
% timesteps from the ode solver.
% later we need to make sure that bubble graph-ing still works, which is in the
% manual integration loop originally

% this code builds upon my "sgLV_v03e2.m" implementation by enabling
% easy maintenance of Target species seeding while varying their initial cell
% densities and maintenance of Target species initial global cell density
% while varying their seeding, as needed

% ON 09/21/2023, I believe I have made the code scalable to any number of
% target species, except for when we need to keep the initial target cell
% densities the same but vary their seeding positions. I WILL HAVE TO FIX
% AT A LATER TIME

% This coe builds on "sgLV_v03f.m" by enabling conversion of a previosly
% considered background species to now consideration of a target species,
% when 'newTarget = 1'

% currently everyrhing works as intended. But, newly converted target
% species will have about 40-times less grid positions than the originally
% selected two target species.

% This code builds on "sgLV_v03 series" by considering all species
% background and then converting them to target species as needed. This way
% a grid position can be occupied only by a background species or a target
% species:
% it worked when I tried for 2 target species.
% next, confirm that it works for greater number of target species too,
% before delving into many many simulations of it.

%%% 16-Feb-2024 %%%
% it worked for 3 and 5 target species too. Next I will test for 9 target species.


clear
close all
clc

%%%% FIRST CHAPTER: RUN ONLY THE FIRST CHAPTER UNTIL Line ~200 WHEN
%%%% CHANGING A NEW # OF TARGET SPECIES

set(0,'DefaultFigureVisible','off')

taskID = str2double(getenv('SLURM_ARRAY_TASK_ID'));

rng(taskID)
%rng('shuffle')

day = yyyymmdd(datetime);

mkdir(num2str(day))
mkdir([num2str(day) '/Global'])

Ntot = 100;

% set to "1" and "run = 0" and "taskID = 0" for new seeding and initialization of the species
init = 0;
run = 0;

% if init = 1, setting this to "1" will load a previous sedding and initialization
% while enabling to change the selected number of target species. Alternatively, setting this
% to "0" will mean a brand-new community seeding and initialization
changeNf = 1;

% if "init = 1" is run once, the initialization folder should be updated
% before the rest of the code can correctly handled.

% current # of target species
Nf = 9;

% current # of background species
NNat = Ntot - Nf;

newrun = taskID * 1;

% filename = [pwd '/20230821/run_' num2str(run)];
filename = [pwd '/initialization/run_' num2str(run)];
filename_new = [pwd '/' num2str(day) '/run_' num2str(newrun)];
filename_new_global = [pwd '/' num2str(day) '/Global/run_' num2str(newrun)];

Cmin = 0; %1.0e-9;
    
n1 = 30;  % Number of rows
n2 = 30;  % Number of columns
    
% community size is n1-by-n2
n = n1 * n2;


if init == 1
%%%% first initialization

    if changeNf == 0

        mean_val_Gn = 0.6; % 0.6 is default
        std_val_Gn = 0.3; % 0.3 is default
        
        % Assign and fix seeding cell densities for the Background species
        CiniNat = rand(n1, n2);
            l6 = length(find(CiniNat <= Cmin));
            while l6 > 0
                CiniNat(CiniNat <= Cmin) = rand(1, l6); %0
                l6 = length(find(CiniNat <= Cmin));
            end
        Cini = CiniNat * 1;
        save([filename_new '_Cini.mat'])
    
        % Generate random species assignments
        Nat_speciesAssignments = randi(Ntot, [n1 n2]);  % Randomly assign species to grid sites
        speciesAssignments = Nat_speciesAssignments * 1;
        save([filename_new '_speciesAssignments.mat'])
    
        muNatVec = normrnd(mean_val_Gn, std_val_Gn, Ntot, 1); 
            l2 = length(find(muNatVec < 0));
            while l2 > 0
                muNatVec(muNatVec < 0) = normrnd(mean_val_Gn, std_val_Gn, l2, 1); %0;
                l2 = length(find(muNatVec < 0));
            end
    
        muVec = muNatVec * 1;
        save([filename_new '_muVec.mat'])
    
    elseif changeNf == 1

        load([filename '_Cini.mat'], 'Cini')
        load([filename '_speciesAssignments.mat'], 'speciesAssignments')
        load([filename '_muVec.mat'], 'muVec')

    end
        
%     muNfVec = [muNatVec(end-(Nf-2)+1:end); muNfVec];
%     muNatVec = muNatVec(1:end-(Nf-2));

%     muNfVec = [muVec(end-(Nf-2)+1:end); muVec];
%     muNatVec = muVec(1:end-(Nf-2));

    muNfVec = muVec(end-(Nf)+1:end);
    muNatVec = muVec(1:end-(Nf));
    
    Nfno = n / Nf;
    % create a pre-assignment for list for target species in the proportion
    % Nfno defined in the line above. Later this list will be randomly shuffled for
    % each random assignment attempt
   
    Nflist = zeros(1, Nf * Nfno);

    for i = 1:Nf

        Nflist(1, Nfno*(i-1)+1 : Nfno * i) = ones(1, Nfno)*(NNat + i);

    end

    Nf_speciesAssignments = reshape(Nflist, [n1, n2]);
    
    CiniNf = zeros(n1, n2);
    
    for i = 1:Nf
    
%         [rf, cf] = find(Nat_speciesAssignments == (NNat + i));
        [rf, cf] = find(speciesAssignments == (NNat + i));

        for j = 1:length(rf)
    
%             Nf_speciesAssignments(rf(j), cf(j)) = Nat_speciesAssignments(rf(j), cf(j));
            Nf_speciesAssignments(rf(j), cf(j)) = speciesAssignments(rf(j), cf(j));
%             CiniNf(rf(j), cf(j)) = CiniNat(rf(j), cf(j));
            CiniNf(rf(j), cf(j)) = Cini(rf(j), cf(j));
%             CiniNat(rf(j), cf(j)) = 0;
            Cini(rf(j), cf(j)) = 0;
        end
    
    
    end
    
    Nat_speciesAssignments = speciesAssignments * 1;
    CiniNat = Cini * 1;

    save([filename_new '_' num2str(Nf) '_CiniNat.mat'])
    save([filename_new '_' num2str(Nf) '_CiniNf.mat'])
    save([filename_new '_' num2str(Nf) '_Nat_speciesAssignments.mat'])
    save([filename_new '_' num2str(Nf) '_Nf_speciesAssignments.mat'])
    save([filename_new '_' num2str(Nf) '_muNatVec.mat'])
    save([filename_new '_' num2str(Nf) '_muNfVec.mat'])
    
end

% %%% UNCOMMENT WHEN CHANGING THE NUMBER OF TARGET SPECIES !!! %%%
% fprintf("Is the initialization folder up-to-date with the matrices to be loaded? \n\n")
% fprintf("If so press any button and comment out these lines next time, if not stop and take care of it \n\n")
% pause;



%%% SECOND CHAPTER %%%
% make BG=1 if you want to plot bubble graphs every tint and record a
% timelapse movie of them.
BG = 0; 


tic

% make ss=1 if you want to load a previously saved initialization of the
% Background speices
ss = 1; 

% make sf=1 if you want to load a previously saved initialization of the
% Target species
sf = 1; 
% make varyCiniNf=1 if you want to randomly assign new initial cell
% densities to Target species while maintaining their previous seeding
% locations and still from the same distribution
% make varyCiniNf=2 if you want to randomly assign new initial cell
% densities to Target species from a randomly chosen subrange of the original distribution
% while maintaining their previous seeding locations

% make varyCiniNf=3 if you want to randomly assign new initial cell
% densities to Target species from a global randomly chosen density from U(0, mfi) then distributed by rnd(1, mfi)
% into the local populations while maintaining their previous seeding locations
varyCiniNf = 3;

% make varyNfseeding=1 if you want to randomly assign new seeding positions 
% to Target species while maintaining their previous global cell density
% make varyNfseeding=2 if you want to randomly assign new seeding positions
% to Target species by using the shuffled version of a previously used
% initial cell density list
varyNfseeding = 0;

% make sInt=1 if you want to load a previously saved interaction matrix
sInt = 1; 
mean_val_intr = 0; % 0 is default
std_val_intr = 0.5; % 0.5 is default

% make sCs=1 if you want to load previously used site colormap
sCs = 1; 

% make sCsp=1 if you want to load previously used species colormap
sCsp = 1; 

% make sNatmu=1 if you want to load a previously created Background species growth rates vector
sNatmu = 1;

% make sNfmu=1 if you want to load a previously used growth rates of the Target species
sNfmu = 1;
mean_val_Gf = 0.6; % 0.6 is default
std_val_Gf = 0.3; % 0.3 is default


if sCs == 0
    % create a color map with unique color to represent the distinct
    % Background or Target species sittig at each grid site
    site_colormap = makeSiteColorMap(n1, n2); % first n1*n2 # rows are for the Background species of the sites and the each following n1*n2 # rows are for the Target species of the sites
    save([filename_new '_site_colormap.mat'])

elseif sCs == 1

    load([filename '_site_colormap.mat'], 'site_colormap')
    site_colormap = site_colormap(:, 1:3); % to ensure compatibility with Python generated colormaps

end

if sCsp == 0
    % create a color map with unique color to represent each distinct
    % species 
    species_colormap = makeSpeciesColorMap(Ntot); % first NNat # rows are for the Background species and the rest is for the Target species
    save([filename_new '_species_colormap.mat'])

elseif sCsp == 1

    load([filename '_species_colormap.mat'], 'species_colormap')
    species_colormap = species_colormap(:, 1:3); % to ensure compatibility with Python generated colormaps

end

if sNfmu == 0
    % create a target species growth rates vector 
    muNfVec = normrnd(mean_val_Gf, std_val_Gf, Nf, 1);
    l1 = length(find(muNfVec < 0));
    while l1 > 0
        muNfVec(muNfVec < 0) = normrnd(mean_val_Gf, std_val_Gf, l1, 1); %0;
        l1 = length(find(muNfVec < 0));
    end
    save([filename_new '_muNfVec.mat'])


elseif sNfmu == 1

    load([filename '_' num2str(Nf) '_muNfVec.mat'], 'muNfVec')

end

if sNatmu == 0

    muNatVec = normrnd(mean_val_Gn, std_val_Gn, NNat, 1); 
    l2 = length(find(muNatVec < 0));
    while l2 > 0
        muNatVec(muNatVec < 0) = normrnd(mean_val_Gn, std_val_Gn, l2, 1); %0;
        l2 = length(find(muNatVec < 0));
    end
    save([filename_new '_muNatVec.mat'])

elseif sNatmu == 1

    load([filename '_' num2str(Nf) '_muNatVec.mat'], 'muNatVec')

end


% scale the nearest neighbor distance
dmin = 1;

gridDistances = calcEucDist(n1, n2, dmin);

% death rate
delta = 0.05;

% define simulation time interval
tend = 50;
tspan = [0 tend];
tint = 5; %tend; %10;


if ss == 0

    % Assign and fix seeding cell densities for the Background species
    CiniNat = rand(n1, n2);
%     l6 = length(find(CiniNat <= Cmin));
%     while l6 > 0
%         CiniNat(CiniNat <= Cmin) = rand(1, l6); %0
%         l6 = length(find(CiniNat <= Cmin));
%     end
    save([filename_new '_CiniNat.mat'])

    % Generate random species assignments
    Nat_speciesAssignments = randi(NNat, n1, n2);  % Randomly assign species to grid sites
    save([filename_new '_Nat_speciesAssignments.mat'])

elseif ss == 1

    load([filename '_' num2str(Nf) '_CiniNat.mat'], 'CiniNat')
    load([filename '_' num2str(Nf) '_Nat_speciesAssignments.mat'], 'Nat_speciesAssignments')

end

if sf == 0

    frc = 0.5;

    mfi = floor(frc * n1 * n2); % number of grid sites seeded with Target species "i"

    % Assign and fix seeding cell densities for the target species
    CiniNf = zeros(n1, n2);
    
    % Generate random indices for replacement
    indices = randperm(n1 * n2, mfi);
    row_f = mod(indices - 1, n1) + 1;
    col_f = ceil(indices / n1);
    
    % Generate random numbers from uniform distribution
    replacement_values = rand(1, mfi);
    l3 = length(find(replacement_values <= Cmin));
    while l3 > 0
        replacement_values(replacement_values <= Cmin) = rand(1, l3); %0
        l3 = length(find(replacement_values <= Cmin));
    end
    
    % Replace selected elements with random values
    CiniNf(indices) = replacement_values; 
    save([filename_new '_CiniNf.mat'])

    shuffledNfList = Nflist(randperm(numel(Nflist)));
    Nf_speciesAssignments = zeros(n1, n2);
    Nf_speciesAssignments(:) = shuffledNfList;
%     Nf_speciesAssignments = NNat + randi(Nf, n1, n2);  % Randomly assign species to grid sites
    save([filename_new '_Nf_speciesAssignments.mat'])

elseif sf == 1

    load([filename '_' num2str(Nf) '_CiniNf.mat'], 'CiniNf')
    load([filename '_' num2str(Nf) '_Nf_speciesAssignments.mat'], 'Nf_speciesAssignments')

    %%%%% previously non-scalable code snippet; test if it is scalable and
    %%%%% working now
    for ti = 1:Nf  

            f_globalCap = 0; % use this below to store the current total biomass of a target species ti

            [vf, idx_f] =  find(Nf_speciesAssignments == NNat + ti);
            f_frozen = zeros(1, length(idx_f)) ;  % use this below to store in a vector the current local cell densities for a target species ti

                for j=1:length(idx_f)

                    f_frozen(j) = CiniNf(vf(j), idx_f(j));
                    f_globalCap = f_globalCap + CiniNf(vf(j), idx_f(j));
    
                end

    end
    %%%%%

%%%%% need to modify according to the new initialization scheme
%     if varyCiniNf == 1 % randomly change the initial cell densities for the target species by keeping their seeding the same
%        
%         rnd_replace = rand(1, mfi);
%         l4 = length(find(rnd_replace <= Cmin));
%         while l4 > 0
%             rnd_replace(rnd_replace <= Cmin) = rand(1, l4); %0
%             l4 = length(find(rnd_replace <= Cmin));
%         end
%         CiniNf(CiniNf > 0) = rnd_replace;
%         save([filename_new '_CiniNf.mat'])
% 
%     end
        
%%%%% need to modify according to the new initialization scheme
%     if varyCiniNf == 2 % this also does more or less the same as 'varyCiniNf' above
%            
%         a_new = rand();
%         b_new = rand();
%         
%         if a_new > b_new
%            
%             ap = a_new;
%             a_new = b_new;
%             b_new = ap;
%             clear ap
%         
%         end
%   
%         rnd_replace = a_new + (b_new - a_new) * rand(1, mfi);
% 
%         l42 = length(find(rnd_replace <= Cmin));
%         while l42 > 0
%             rnd_replace(rnd_replace <= Cmin) = rand(1, l42); %0
%             l42 = length(find(rnd_replace <= Cmin));
%         end
%         CiniNf(CiniNf > 0) = rnd_replace;
%         save([filename_new '_CiniNf.mat'])
% 
%     end

    %%%%% previously non-scalable code snippet; test if it is scalable and
    %%%%% working now

    if varyCiniNf == 3 % this assigns new initial cell densities to the same seeding positions by sampling a total (global) cell density
% for a given target species and randomly distributing it into the local
% populations
        Nf_speciesAssignmentsVec = Nf_speciesAssignments(:);
        CiniNfvec = CiniNf(:);

        for ti = 1:Nf

            rng('shuffle')

            linx_sf = find(Nf_speciesAssignmentsVec == NNat + ti);

            [vsf_, ix_sf] =  find(Nf_speciesAssignments == NNat + ti);

            a_new_s = 0;
            b_new_s = length(nonzeros(CiniNfvec(linx_sf)));

            s = a_new_s + (b_new_s - a_new_s) * rand();

            d_s = rand(1, b_new_s);

            s_vec = ( d_s / sum(d_s) ) * s;

            ks = 1;

            for j=1:length(ix_sf)
                
                if CiniNf(vsf_(j), ix_sf(j)) > 0

                    CiniNf(vsf_(j), ix_sf(j)) = s_vec(ks);

                    ks = ks + 1;
        
                end

            end
        

        end
            
%             save([filename_new '_CiniNf.mat'])
    
    end

    %%%%%

%     %%%%% BROKEN PIECE OF CODE
%         %%% f_globalCap is the global sum of the original target
%         %%% species initial cell densities. As the code was attempted to
%         %%% made scalable the indices of f1_globalCap and f2_globalCap have
%         %%% been removed. So the portion below will not work and need a
%         %%% fix. Hence I commented them out.
%     if varyNfseeding == 1 % randomly vary both seeding poisitons and initial cell densities
%         shuffledNfList = Nflist(randperm(numel(Nflist)));
%         Nf_speciesAssignments = zeros(n1, n2);
%         Nf_speciesAssignments(:) = shuffledNfList;
%         save([filename_new '_Nf_speciesAssignments.mat'])
% 
%          % Generate random indices for replacement
%         indices = randperm(n1 * n2, mfi);
%         
%         % Generate random numbers from normal distribution
%         rnd_replace2 = rand(1, mfi);
%         l5 = length(find(rnd_replace2 <= Cmin));
%         while l5 > 0
%             rnd_replace2(rnd_replace2 <= Cmin) = rand(1, l5); %0
%             l5 = length(find(rnd_replace2 <= Cmin));
%         end
% 
% 
%         % Assign and fix seeding cell densities for the Target species
%         CiniNf = zeros(n1, n2);
%         CiniNf(indices) = rnd_replace2; 
%     
%         %%%%% previously non-scalable code snippet; test if it is scalable and
%         %%%%% working now
%         
%         %%% f_globalCap is the global sum of the original target
%         %%% species initial cell densities. As the code was attempted to
%         %%% made scalable the indices of f1_globalCap and f2_globalCap have
%         %%% been removed. So the portion below will not work and need a
%         %%% fix. Hence I commented them out.
%         for ti = 1:Nf
%             
%             foc_global = 0;
%             [vf_, ix_f] =  find(Nf_speciesAssignments == NNat + ti);
% 
%             for j=1:length(ix_f)
% 
%                 foc_global = foc_global + CiniNf(vf_(j), ix_f(j));
%     
%             end
% 
%             for j=1:length(ix_f)
% 
%                 CiniNf(vf_(j), ix_f(j)) = CiniNf(vf_(j), ix_f(j)) / (foc_global / f_globalCap);
%         
%             end
% 
% 
%         end
%         %%%%%
% 
%         save([filename_new '_CiniNf.mat'])
% 
%     end


%      if varyNfseeding == 2
%         
%         shuffledNfList = Nflist(randperm(numel(Nflist)));
%         Nf_speciesAssignments = zeros(n1, n2);
%         Nf_speciesAssignments(:) = shuffledNfList;
%         save([filename_new '_Nf_speciesAssignments.mat'])
% 
%         %%%%% non-scalable code snippet
%         shuffledf1_frozen = f1_frozen(randperm(numel(f1_frozen)));
%         shuffledf2_frozen = f2_frozen(randperm(numel(f2_frozen)));
% 
%         [vf1_, ix_f1] =  find(Nf_speciesAssignments == NNat + 1);
%         [vf2_, ix_f2] =  find(Nf_speciesAssignments == NNat + 2);
%         %%%%%
% 
%          % Generate random indices for replacement
%         indices = randperm(n1 * n2, mfi);
%         
%         % Generate random numbers from a uniform distribution
%         rnd_replace2 = rand(1, mfi);
%         l5 = length(find(rnd_replace2 <= Cmin));
%         while l5 > 0
%             rnd_replace2(rnd_replace2 <= Cmin) = rand(1, l5); %0
%             l5 = length(find(rnd_replace2 <= Cmin));
%         end
% 
%         % Assign and fix seeding cell densities for the Background species
%         CiniNf = zeros(n1, n2);
%         CiniNf(indices) = rnd_replace2; 
% 
%         %%%%% non-scalable code snippet
%         for j=1:length(ix_f1)
% 
%             CiniNf(vf1_(j), ix_f1(j)) = shuffledf1_frozen(j);
%     
%         end
%     
%         for j=1:length(ix_f2)
% 
%             CiniNf(vf2_(j), ix_f2(j)) = shuffledf2_frozen(j);
%     
%         end
%         %%%%%%
% 
%         save([filename_new '_CiniNf.mat'])
% 
%     end
%%%%%%%%

end

% if ~isequal(mfi, length(find(CiniNf > 0)))
%     fprintf('number of sites inoculated with Target species is not as intended \n\n')
% end


% maximum specific growth rates
muVec = [muNatVec(:); muNfVec(:)];

if sInt == 0

    % create the interaction matrices by complying with Zhengqing' s initialization
  % that is between any two species interactions are single-edged and can be either positive or negative but not both      
      gammaInteraction = normrnd(mean_val_intr, std_val_intr, [Ntot, Ntot]);

      % assuming the rows 1 to NNat and columns 1 to NNat represent the
      % Background species and the rest represent the Target species
        
        a = gammaInteraction; q = a;
        
        a(gammaInteraction < 0) = 0; % positive interactions
        q(gammaInteraction > 0) = 0; q = -1*q; % negative interactions
        
        %%% to comply with Zhengqing' s implementation make the
        %%% interactions single-edged
         [a, q] = makeIntSingleEdged(a, q);

         % Combined interaction matrix
        interactionMatrix = a - q;  % Replace this with your actual interaction matrix

        save([filename_new '_gammaInteraction.mat'])
        save([filename_new '_interactionMatrix.mat'])

elseif sInt == 1

        load([filename '_gammaInteraction.mat'], 'gammaInteraction')
        load([filename '_interactionMatrix.mat'], 'interactionMatrix')

end


[interactionVectors, Nat_growthRateMatrix, Nf_growthRateMatrix, gridDistances] = CalcInteractions(n1, n2, gridDistances, Nat_speciesAssignments, Nf_speciesAssignments, interactionMatrix, muVec);

% set up the interaction matrices after scaling by the position
[F, G] = formMatrix04(n1, n2, interactionVectors);

% initialize
growthRateVec = [Nat_growthRateMatrix(:); Nf_growthRateMatrix(:)];
GR = growthRateVec * 1;

% initial cell density state vector
y0 = [CiniNat(:); CiniNf(:)];
        
t = 0;

if BG == 1
    
    % Create a VideoWriter object
    video_filename = [filename_new '_bubble_graph_video.mp4'];
    frame_rate = 2;  % Set the desired frame rate (frames per second)
    video = VideoWriter(video_filename, 'MPEG-4');
    video.FrameRate = frame_rate;  % Set the frame rate
    open(video);

    sz = floor(5e2/(n1*n2/36));
    lw = 1;

    V1 = y0(1:n1*n2, 1);
    V2 = y0(n1*n2+1:end, 1);

    V1(V1 == 0) = nan;
    V2(V2 == 0) = nan;

    % Create bubble graphs for this state
    fig1 = figure('Color', 'white'); % Set the figure background color to white
    subplot(1,2,1)
    BubbleGraph(n1, n2, V1, Nat_speciesAssignments, site_colormap, species_colormap, sz, lw);
    title('Background species', 'FontSize', 22)

    subplot(1,2,2)
    BubbleGraph(n1, n2, V2, Nf_speciesAssignments, site_colormap, species_colormap, sz, lw);
    title('Target species', 'FontSize', 22)

    sgtitle(['t = ' num2str(t)], 'FontSize', 36);
    drawnow;

    % Capture the current figure as a frame and write to the video
    frame = getframe(gcf);
    writeVideo(video, frame);

    % Clear the figure for the next frame
    clf;     

end

% solve using a built-in numerical integrator
   
    opts = odeset('NonNegative', 1);       
    sol = ode15s(@(t, y)FUNsgLV(t, y, delta, F, G, GR), tspan, y0, opts);
    time = 0:1:tend;
    y = deval(sol, time);
    cellDens = y' * 1;
    TimeTrace = [time' cellDens];
    pidx = 1;

%         % hand-written Euler solver dt = 0.1;
%         
%         time = zeros(tend/dt + 1, 1); cellDens = zeros(tend/dt + 1,
%         2*n1*n2);
%         
%         time(1, :) = t; cellDens(1, :) = y0;
%         
%          y = y0;
%         
%            for i = 1 : tend/dt + 1
%         
%                dy_dt = GR .* y .* (1 - y) - (delta .* y ./ (F*y + 1) ) -
%                (G*y .* y) ;
%         
%                y = y + dy_dt * dt; y(y < Cmin) = 0; t = i*dt;
%                     
%                time(1+i, :) = t; cellDens(1+i, :) = y;
%                 
%                pidx = 1;
%                
%                if (BG == 1) && (mod(i*dt, tint) == 0)
% 
%                         V1 = y(1:n1*n2, 1); V2 = y(n1*n2+1:end, 1);
%                                            
%                         V1(V1 == 0) = nan; V2(V2 == 0) = nan;
%                                 
%                         % Create bubble graphs for this state fig2 =
%                         figure('Color', 'white'); % Set the figure
%                         background color to white subplot(1,2,1)
%                         BubbleGraph(n1, n2, V1, Nat_speciesAssignments,
%                         site_colormap, species_colormap, sz, lw);
%                         title('Background species', 'FontSize', 22)
%                        
%                         subplot(1,2,2) BubbleGraph(n1, n2, V2,
%                         Nf_speciesAssignments, site_colormap,
%                         species_colormap, sz, lw); title('Target species',
%                         'FontSize', 22) sgtitle(['t = ' num2str(t)],
%                         'FontSize', 36); drawnow;
%                       % Capture the current figure as a frame and write
%                       to the video
%                          frame = getframe(gcf); writeVideo(video, frame);
%     
%                          % Clear the figure for the next frame clf;
%         %             pause;
%                end
%                 pidx = pidx + 1;
%            end
%  
%            if BG == 1
%               % Close the video
%                 close(video);
%            end
%            
%         close all
% 
%          TimeTrace = [time cellDens];
            
fprintf(['Simulation (not including visualization or writing output .csv files) took ' num2str(toc/60) ' min \n\n'])
               
           fig3 = figure(pidx);
           subplot(2,1,1)
                for site_index = 1:n1*n2
                    site_color = site_colormap(site_index, :);
                    plot(time, cellDens(:, site_index), 'Color', site_color, 'LineWidth', 3);
                    title('Background species local cell density')
%                     xlabel('Time');
%                     ylabel('Local Cell Density');
                    set(gca, 'FontSize', 18);
                    hold on;
                end
                hold off

           subplot(2,1,2)
                for site_index = 1:n1*n2
                    site_color = site_colormap(site_index, :);
                    plot(time, cellDens(:, n1 * n2 + site_index), 'Color', site_color, 'LineWidth', 3);
                    title('Target species local cell density')
                    xlabel('Time');
%                     ylabel('Local Cell Density');
                    set(gca, 'FontSize', 18);
                    hold on;
                end
                hold off
        
           sgtitle({['Square lattice, n_1 = ' num2str(n1) ', n_2 = ' num2str(n2) ], ['d_{min} = ' num2str(dmin) ', \delta = ' num2str(delta)]}, 'FontSize', 22)
           drawnow;
           saveas(fig3, [filename_new '_LocalCellDensities.pdf'])

           writematrix(TimeTrace, [filename_new '_LocalCellDensitiesTimeCourse.csv'])
  

globalNat = zeros(length(time), NNat);
globalNf = zeros(length(time), Nf);

Nat_speciesAssignmentsVec = Nat_speciesAssignments(:);
Nf_speciesAssignmentsVec = Nf_speciesAssignments(:);

 for i = 1:NNat

    Natix = find(Nat_speciesAssignmentsVec == i);

    sumNat = 0;

    for j= 1:length(Natix)

        sumNat = sumNat + cellDens(:, Natix(j));

    end

        globalNat(:, i) = sumNat;


 end


 for i = 1:Nf

    Nfix = find(Nf_speciesAssignmentsVec == (NNat + i));

    sumNf = 0;

    for h= 1:length(Nfix)

        sumNf = sumNf + cellDens(:, n1 * n2 + Nfix(h));

    end

        globalNf(:, i) = sumNf;
        

 end

           
          fig4 = figure(pidx + 1);
          subplot(2,1,1)
                for Natspecies_index = 1:NNat
                    species_color = species_colormap(Natspecies_index, :);
                    plot(time, globalNat(:, Natspecies_index), 'Color', species_color, 'LineWidth', 3);
                    title('Background species global cell density')
%                     xlabel('Time');
%                     ylabel('Local Cell Density');
                    set(gca, 'FontSize', 18);
                    hold on;
                end
                hold off

           subplot(2,1,2)
                for Nfspecies_index = 1:Nf
                    species_color = species_colormap(NNat + Nfspecies_index, :);
                    plot(time, globalNf(:, Nfspecies_index), 'Color', species_color, 'LineWidth', 3);
                    title('Target species global cell density')
                    xlabel('Time');
%                     ylabel('Local Cell Density');
                    set(gca, 'FontSize', 18);
                    hold on;
                end
                hold off
        
           sgtitle({['Square lattice, n_1 = ' num2str(n1) ', n_2 = ' num2str(n2) ], ['d_{min} = ' num2str(dmin) ', \delta = ' num2str(delta)]}, 'FontSize', 22)
           drawnow;
           saveas(fig4, [filename_new_global '_GlobalCellDensities.pdf'])


OutputTable1 = [time' globalNat globalNf];
writematrix(OutputTable1, [filename_new_global '_GlobalCellDensitiesTimeCourse.csv'])