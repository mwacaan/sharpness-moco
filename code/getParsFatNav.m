pars.region = [-3, 0, 3, 4]; % [min_R1 max_R1 max_R2 max_R3]. Note: defined with coordinate offset (legacy)
pars.smf = 3; % smoothing factor for Gaussian blur
pars.mrelcb = .5; % maximum size of relative confidence bounds for valid cluster
pars.wI = 1; % pixel intensity weights for clustering (relative to coordinates)
pars.numpoints = 500; 
pars.numIcl = 4; % number of intensity clusters
pars.border_distance = 1.1;
