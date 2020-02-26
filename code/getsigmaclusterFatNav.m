function [ better_signed, worse_signed, brd_crds, brd_ind, AddRes ] = getsigmaclusterFatNav( map_corr, data_corr, data_uncorr, pars, field_name, plotsavedir )
% Calculate sigma for every point near border.
warning off 

mapInf=1e3; % map padding zeros saturation value

%% Preparation for clustering
% Load parameters
border_distance = pars.border_distance;
min_R1 = pars.region(1);
max_R1 = pars.region(2);
max_R2 = pars.region(3);
max_R3 = pars.region(4);
k_R3 = pars.numIcl;
cluster_size = pars.numpoints;

mask=and(data_corr>0,data_uncorr>0);
frac=sum(mask(:))/numel(mask);
disp(['Union mask of data_corr and data_uncorr percentage: ' num2str(frac)])
% mask data
data_corr(~mask)=0;
data_uncorr(~mask)=0;

% Load and scale data. Clustering is done on corrected data. Can be changed
% to uncorrected data.
map = map_corr;
data = data_corr;
map(~mask)=mapInf; %saturate padding zeroes in map

% Define scaling factors and scale data
scale_I = pars.wI./prctile(nonzeros(data.*(map<max_R3)),90);
scale_c = max(size(map)); % used later on
data = scale_I.*data;

% Define regions with distance transform values. The map variables are
% binary.
R1_map = map>min_R1 & map<=max_R1;
R1_val = data.*R1_map;
R2_map = map>max_R1 & map<max_R2;
R2_val = data.*R2_map;
R3_map = map>max_R2 & map<max_R3;
R3_val = data.*R3_map;
% DIPimage smoothing
%R3_val_sm = double(gaussf(data,pars.smf)).*R3_map;
% matlab smoothing
R3_val_sm = double(imgaussian(data,pars.smf)).*R3_map;

I_R3_sm = R3_val_sm(R3_val~=0);

% Get (scaled) spatial coordinates of ROI regions
[x_R1, y_R1, z_R1] = ind2sub(size(map),find(R1_val~=0));
coords_R1 = cat(2, x_R1, y_R1, z_R1)./scale_c;
[x_R2, y_R2, z_R2] = ind2sub(size(map),find(R2_val~=0));
coords_R2 = cat(2, x_R2, y_R2, z_R2)./scale_c;
[x_R3, y_R3, z_R3] = ind2sub(size(map),find(R3_val~=0));
coords_R3 = cat(2, x_R3, y_R3, z_R3)./scale_c;

%% Clustering
% First kmeans-clustering: on intensity value
ind_R3 = kmeans(cat(2,coords_R3,I_R3_sm),k_R3, 'MaxIter', 1000, 'Replicates', 10);
ind_R3subcl = ind_R3; % used for subclustering. 

% Second k-means clustering: cluster intensity clusters on spatial
% coordinates. The number of subclusters is determined by the clustersize.
k_subcl = zeros(1,k_R3);
C_reps = [];
I_inds = [];

for i = 1:k_R3
    numpoints = length(nonzeros(ind_R3==i));
    k_subcl(i) = ceil(numpoints./cluster_size);
    [ind_rep,C_rep] = kmeans(coords_R3(ind_R3==i,:),k_subcl(i), 'MaxIter', 1000, 'Replicates', 10);
    
    % Correct the index values to avoid subclusters of different intensity
    % clusters having the same index.
    valcor = sum(k_subcl(1:i))-k_subcl(i);
    ind_rep = ind_rep+valcor;
    ind_R3subcl(ind_R3==i)=ind_rep;
    
    C_reps = cat(1,C_reps,C_rep);
    I_inds = cat(2,I_inds,i.*ones(1,k_subcl(i)));
end
k_tot = sum(k_subcl);
tempind = 1:k_tot; % when ROI is very small, use tempind = I_inds

% Propagate clustering to R1 and R2, get border coordinates and indices.

R2toR3 = knnsearch(coords_R3,coords_R2);
ind_R2 = ind_R3subcl(R2toR3);
R1toR2 = knnsearch(coords_R2,coords_R1);
ind_R1 = ind_R2(R1toR2);
border = double(dt(map<0)<border_distance & dt(map<0)>0);
brd_crds = find(border~=0)';
indsR1 = find(R1_val~=0);
brd_ind = ind_R1(ismember(indsR1,brd_crds));

%% Error function fitting
% Use loop over corrected and uncorrected data to do fitting for both sets.
% Same clustering is used for both sets of data.

sigma_comp = cell(1,2);
err_comp = cell(1,2);
valid_comp = cell(1,2);
n_invalid = zeros(1,2);
aggregates = cell(1,2);
for p = 1:2
    map = map_corr; % map_corr is used here, as it was used for clustering as well. map_corr and map_uncorr are expected to be nearly identical.
    if p==1
        data = data_corr;
        ctag='corr'
    else
        data = data_uncorr;
        ctag='uncorr'
    end
    map(~mask)=mapInf;

    % Scale data, get intensity and distance transform data for the
    % regions.
    data = scale_I.*data;

    R1_map = map>min_R1 & map<=max_R1;
    I_R1 = data(R1_map);
    DT_R1 = map(R1_map);
    
    R2_map = map>max_R1 & map<max_R2;
    I_R2 = data(R2_map);
    DT_R2 = map(R2_map);
    
    R3_map = map>max_R2 & map<max_R3;
    I_R3 = data(R3_map);
    DT_R3 = map(R3_map);

    % Create profiles used for fitting
    prof_fit = cell(1,k_tot);
    for i = 1:k_tot
        try
          Iprof = cat(1,I_R1(ind_R1==tempind(i)),I_R2(ind_R2==i),I_R3(ind_R3subcl==i));
          DTprof = cat(1,DT_R1(ind_R1==tempind(i)),DT_R2(ind_R2==i),DT_R3(ind_R3subcl==i));
          prof_fit{i} = cat(2,DTprof,Iprof);
        catch
          keyboard
        end
    end

    % Do erf fitting
    min_datapoints = cluster_size/10; % 10 is arbitrary. It shouldn't be too high or too low. 
    sigma = zeros(1,k_tot);
    err = zeros(1,k_tot);
    valid = zeros(1,k_tot);
    aggregate = zeros(1,k_tot);
    f_fit = cell(1,k_tot);
    for i = 1:k_tot  
        [f_fit{i}, sigma(i), cb, diff, ~, sse] = errorfunctionFATNAV(prof_fit{i});
        l_R1 = length(nonzeros(ind_R1==i));
        err(i) = abs(cb(1)-sigma(i));
        cb_rel = 2*err(i)./sigma(i);
        sse_rel = sse/l_R1;
        
        % validity testing based on 3 criteria, build binary coding for
        % invalidity
        c1=abs(diff)>2*sse_rel; 
        c2=abs(cb_rel)<pars.mrelcb;
        c3=l_R1>min_datapoints;
        if c1 && c2 && c3
            valid(i) = 1;
            aggregate(i)=0;
        else
            aggregate(i)=~c1*-1 + ~c2*-2 + ~c3*-4;
        end
    end

    % This plots the first 20 profiles and fit lines, which 
    % should be sufficient to see whether something is wrong. 
    colors=[0 0 0;0 0 1];
    lcolors=[1 0 0;0 0 0];
    f = figure('visible', 'off');
    if k_tot>19
        i_max = 20;
    else
        i_max = k_tot;
    end
    for i = 1:i_max
        subplot_tight(4,5,i,.02)
        plot(prof_fit{i}(:,1),prof_fit{i}(:,2),'.','color',colors(valid(i)+1,:), ...
          'markersize',2)
        hold on
        h=plot(f_fit{i});
        set(h,'color',lcolors(valid(i)+1,:))
        axis([-3 4 -.5 1.5])
        legend('hide')
        xlabel(''); ylabel('')
    end
    name = [field_name, '_',ctag,'_profiles'];
    filename = fullfile(plotsavedir,[name '.png']);
    disp(['Creating ' filename])
    print(filename,'-dpng','-r300')
    close(f)
  
    % Save variables to compare corrected and uncorrected data
    sigma_comp{p} = sigma;
    err_comp{p} = err;
    valid_comp{p} = valid;
    n_invalid(p) = k_tot-length(nonzeros(valid));
    aggregates{p} = aggregate;
end

%% Compare corrected and uncorrected data

validboth = valid_comp{1}.*valid_comp{2};
sc = sigma_comp{1};
ec = err_comp{1};
su = sigma_comp{2};
eu = err_comp{2};
sigma_diff_abs = validboth.*sqrt(abs(su.^2-sc.^2));
signed_indication = 2.*((su.^2-sc.^2)>0)-1; % 1 if better, -1 if worse
sigma_diff_signed = signed_indication.*sigma_diff_abs;

better_signed = validboth.*sigma_diff_signed.*((su-eu)>(sc+ec));
worse_signed = validboth.*sigma_diff_signed.*((sc-ec)>(su+eu));
unc = sqrt(abs(1./(su.^2-sc.^2).*(su.^2.*eu.^2+sc.^2.*ec.^2)));

%% Store additional results 

AddRes.('crdindicesR1') = indsR1; % coordinate indices of R1
AddRes.('indR1') = ind_R1; % cluster indication of subclusters
AddRes.('Iinds') = I_inds; % cluster indication of intensity clusters

AddRes.('sigmacorr') = sigma_comp{1}; % sigma for corrected data
AddRes.('errcorr') = err_comp{1}; % fit error for corrected data
AddRes.('sigmauncorr') = sigma_comp{2}; % sigma for uncorrected data
AddRes.('erruncorr') = err_comp{2}; % fit error for uncorrected data
AddRes.('sigmadiff') = sigma_diff_abs; % absolute squared difference in sigma
AddRes.('sigmadiff_signed') = sigma_diff_signed; % signed squared difference in sigma
AddRes.('uncertainty') = unc; % uncertainty. Not sure if it is relevant

AddRes.('valid') = valid_comp; % 2-cell array of validity indications. cell 1 is corrected, cell 2 uncorrected.
AddRes.('ninvalid') = n_invalid; % number of invalid clusters
AddRes.('validboth') = validboth; % indicates validity of both sets of data
AddRes.('aggregates') = aggregates; % binary indication of invalidity
end

