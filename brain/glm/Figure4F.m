% To run GLM, you will need to request the full dataset from corresponding author Ken Shepard (shepard@ee.columbia.edu). The data files are too large to be included in the Github repository. 

%% Tool paths

% Add nirs-toolbox to path
nirs_toolbox_path = join([getenv('UserProfile') filesep "Documents" filesep "MATLAB" filesep "nirs-toolbox"], "");
addpath(nirs_toolbox_path);
addpath(genpath(nirs_toolbox_path + filesep + "demos"));
addpath(genpath(nirs_toolbox_path + filesep + "external"));

% Add isomesh to path
addpath(join([getenv('UserProfile') filesep 'Documents\MATLAB\iso2mesh'], "")); 

% Custom function
addpath(join([getenv('UserProfile') filesep "Dropbox" filesep "MOANA" filesep "Homer3" filesep "matlab_functions"], ""));


%% Set root directory
root_dir = [getenv('UserProfile') filesep 'Dropbox\MOANA\Python\MOANA3_Python37_Codes\chips\moana3\data\2023-12-04_Petros_Motor6s_MOANA_Pinhole\motor'];

%% load data
raw = nirs.io.loadDirectory([root_dir filesep 'forNT_motor_block_single_source_setting16'], {'subject', 'session', 'run'});

%% Register the probe

% Load the SD file
sd_path = load([getenv('UserProfile') filesep 'Dropbox\MOANA\Homer3\new_moana3_probe\MOANA3_FLEX_4BY4_SOMATOSENSORY_3D.SD'], "-mat");
sd_path = sd_path.SD;

% Grab the anchor list
ancl = sd_path.AnchorList;

% Create coordinate list
Name{1}='';
xyz(1,:)=[0 0 0];
Type{1}='FID-anchor';  % This is an anchor point
Units{1}='mm';

% Now let's add a few more
for i = 1:length(ancl)

    Name{i}=sd_path.AnchorList{i, 2};
    xyz(i,:)=sd_path.DummyPos(i, :, :);
    Type{i}='FID-anchor';  % This is an attractor
    Units{i}='mm';

end

% Create table with these anchors
fid=table(Name',xyz(:,1),xyz(:,2),xyz(:,3),Type',Units',...
    'VariableNames',{'Name','X','Y','Z','Type','Units'});

% Add to probe, register, attach
for i = 1:length(raw)
    probe=raw(i).probe;
    probe.optodes=[probe.optodes; fid];
    probe1020=nirs.util.registerprobe1020(probe);

    lambda=unique(probe1020.link.type);
    fwdBEM=nirs.registration.Colin27.BEM(lambda);
    
    headshape=nirs.registration.getheadshape(fwdBEM.mesh(1));
      
    probe1020=probe1020.register_mesh2probe(fwdBEM.mesh);
    probe1020.defaultdrawfcn='3D mesh';
    
    raw(i).probe = probe1020;
end

%% Remove channels

% Null source-detector separation channels
j = nirs.modules.MOANA_RemoveChannels( );
j.channel_idcs = raw(1).probe.link.source == raw(1).probe.link.detector;
raw_pruned = j.run(raw);

% Source/Detector pairs that are non-functional are set to 0
zero_idcs = true(1, size(raw_pruned(1).data, 2));
for file_idx=1:length(raw_pruned)
    zero_idcs = zero_idcs & (sum(raw_pruned(file_idx).data, 1) == 0);
end
zero_idcs = find(zero_idcs);

% Remove
sources_to_remove = [];
sources_to_inspect = unique(raw_pruned(1).probe.link.source(zero_idcs));
for s=1:length(sources_to_inspect)
    source = sources_to_inspect(s);
    number_of_dead_channels = sum(raw_pruned(1).probe.link.source(zero_idcs) == source);
    if number_of_dead_channels == 30
        sources_to_remove = [sources_to_remove source];
    end
end
j = nirs.modules.MOANA_RemoveChannels( );
j.channel_idcs = find(any(raw_pruned(1).probe.link.source == sources_to_remove, 2));
raw_pruned = j.run(raw_pruned);

% Remove sources and detectors in the middle of the patch, as these do not
% probe relevant depths
j = nirs.modules.MOANA_RemoveChannels( );
j.channel_idcs = find(any(raw_pruned(1).probe.link.source == [6 7 10 11], 2));
raw_pruned = j.run(raw_pruned);

% Remove channels that fall in some intermediate range of distances
j = nirs.modules.MOANA_RemoveChannels( );
j.channel_idcs = find((raw_pruned(1).probe.distances < 23) & (raw_pruned(1).probe.distances > 9));
raw_pruned = j.run(raw_pruned);


%% Rename stims
j = nirs.modules.RemoveStimless( );
j = nirs.modules.RenameStims( j );
j.listOfChanges = {
    'stim_channel1', 'stim'}; %'palm'};
j = nirs.modules.KeepStims( j );
j.listOfStims = {'stim'};
raw_renamed = j.run(raw_pruned);


%% Downsample
j = nirs.modules.Resample( );
Fs = 5;
j.Fs = Fs;
raw_downsampled = j.run(raw_renamed);


%% Trim baseline
j = nirs.modules.TrimBaseline(  );
j.preBaseline  = 30;
j.postBaseline = 30;
raw_trimmed = j.run(raw_downsampled);


%% Convert to optical density and baseline correct
j = nirs.modules.OpticalDensity(  );
od = j.run(raw_trimmed);


%% Convert to hemoglobin.
j = nirs.modules.BeerLambertLaw(  );
j.PPF = [6.3 5.6];
hb = j.run( od );


%% Modify stim durations
stimNames = hb(1).stimulus.keys';
stimNames = stimNames{1};
hb = nirs.design.change_stimulus_duration(hb, stimNames, 6);


%% Create list of trigger and file names
trigger_names = hb(1).stimulus.keys';
file_names = {'001', '002'};


%% Add short separation channels
j = nirs.modules.LabelShortSeperation( );
j.max_distance = 9; 
hb = j.run(hb);
j = nirs.modules.AddShortSeperationRegressors( );
ss_hb = j.run(hb);


%% Subject level pipeline
j = nirs.modules.GLM();
j.verbose = true;

% Baseline will include constant, linear, and second order term
trend_func = @(t) nirs.design.trend.legendre(t, 2);
j.trend_func = trend_func;

% Use short separation channels
j.AddShortSepRegressors = true;

% Specify basis function
j.basis = Dictionary();
b = nirs.design.basis.Canonical();
j.basis('default') = b;

% Run
SubjStats = j.run( hb );

% vizualization
SubjStats(1).probe.defaultdrawfcn = ['2D'];
ob = SubjStats(1).draw('p', [0.2 -0.2]);