%% add fieldtrip to the path and initialize
ftpath = '/home/fsandhaeger/matlab/fieldtrip-20231025/'; % change to your fieldtrip path
addpath(ftpath);
ft_defaults

dsfile = '/home/fsandhaeger/projects/fs_preproc_meg/example_data/meg_example.ds';
bhvfile  = '/home/fsandhaeger/projects/fs_preproc_meg/example_data/bhv_example.mat';

% change to your meg data file
hdr=ft_read_header(dsfile);

%% double check that the recording was continuous and no gaps exist
dt_clock = ft_read_data(dsfile,'chanindx',find(strcmp(hdr.label,'SCLK01')));
cd=diff(dt_clock(:)); % difference between subsequent samples to identify gaps in the recording
cd = cd(cd>0);
cd = abs(((1/hdr.Fs)./cd)-1);
max_dev = max(cd); % maximal deviation from nominal sampling rate

disp(['samples with deviation form nominal sampling rate > 1%: ' num2str(sum(cd>0.01))]);
if max_dev>0.25,
    warning('there are samples with a deviation from the sampling rate of >25%!');
end

%% now all data

% A) load using the fieldtrip preprocessing function
cfg_load = [];
cfg_load.dataset = dsfile;
cfg_load.continuous = 'yes';
data = ft_preprocessing(cfg_load);

% B) if that for some reason doesn't work well, use a more basic function
% dat = ft_read_data(dsfile);
% data = [];
% data.trial{1}  = reshape(dat,size(dat,1),[]); clear dat
% data.label     = hdr.label;
% data.time{1}   = (1:size(data.trial{1},2))/hdr.Fs;
% data.fsample   = hdr.Fs;

idx_meg = find(cellfun(@(x) strcmp(x(1),'M'),data.label));

%% remove channel jumps
%  because this recording doesn't have any channel jumps, we just add some
%  artificially

idx_jump = round(size(data.trial{1},2)*0.34);
data.trial{1}(idx_meg(10),idx_jump:end) = data.trial{1}(idx_meg(10),idx_jump:end)+2*mean(data.trial{1}(idx_meg(10),idx_jump:end));
idx_jump = round(size(data.trial{1},2)*0.89);
data.trial{1}(idx_meg(104),idx_jump:end) = data.trial{1}(idx_meg(104),idx_jump:end)-6*mean(data.trial{1}(idx_meg(104),idx_jump:end));

[data_nojumps jumps] = fs_correct_jumps(data);

figure;
plot(data.trial{1}(idx_meg(10),1:100:end));
hold on;
plot(data_nojumps.trial{1}(idx_meg(10),1:100:end));
title('example channel before and after jump correction');

% resample to 1000Hz
cfg_rs = [];
cfg_rs.resamplefs = 1000;
cfg_rs.detrend = 'no';
data_1000 = ft_resampledata(cfg_rs,data_nojumps);

% bandstop filter to remove line noise
cfg_flt = [];
cfg_flt.bsfilter = 'yes';
cfg_flt.bsfreq = [49 51; 99 101; 149 151];
cfg_flt.bsinstabilityfix = 'split';
data_bs = ft_preprocessing(cfg_flt,data_1000);

% highpass filter
cfg_flt = [];
cfg_flt.hpfilter = 'yes';
cfg_flt.hpfreq = [0.05];
cfg_flt.hpinstabilityfix = 'split';
data_hp = ft_preprocessing(cfg_flt,data_bs);

% lowpass filter
cfg_flt = [];
cfg_flt.lpfilter = 'yes';
cfg_flt.lpfreq = [150];
cfg_flt.lpinstabilityfix = 'split';
data_lp = ft_preprocessing(cfg_flt,data_hp);


% plot some example data
channels_plot = idx_meg(round(linspace(1,length(idx_meg),10))); % select 10 channels
figure;
subplot(2,4,1);
fs_plot_channels(data.time{1}(1:100:end),data.trial{1}(channels_plot,1:100:end));
xlabel('time in s');
ylabel('T');
title('raw data, 10 channels whole exp.');
subplot(2,4,2);
fs_plot_channels(data.time{1}(round(data.fsample*500)+[1:round(data.fsample*60)]),data.trial{1}(channels_plot,round(data.fsample*500)+[1:round(data.fsample*60)]));
xlabel('time in s');
ylabel('T');
title('raw data, 10 channels 60s');
subplot(2,4,5);
dtmp = data.trial{1}(idx_meg,1:10:end);
imagesc(data.time{1}(1:100:end),[],data.trial{1}(idx_meg,1:100:end)-repmat(mean(data.trial{1}(idx_meg,1:100:end),2),1,size(data.trial{1}(:,1:100:end),2)));
caxis([-1 1]*prctile(abs(dtmp(:)),99));
xlabel('time in s');
ylabel('channels');
title('raw data, all channels whole exp.');

subplot(2,4,3);
fs_plot_channels(data_lp.time{1}(1:10:end),data_lp.trial{1}(channels_plot,1:10:end));
xlabel('time in s');
ylabel('T');
title('filtered data, 10 channels whole exp.');
subplot(2,4,4);
fs_plot_channels(data_lp.time{1}(round(data_lp.fsample*500)+[1:round(data_lp.fsample*60)]),data_lp.trial{1}(channels_plot,round(data_lp.fsample*500)+[1:round(data_lp.fsample*60)]));
xlabel('time in s');
ylabel('T');
title('filtered data, 10 channels 60s');
subplot(2,4,7);
dtmp = data_lp.trial{1}(idx_meg,1:10:end);
imagesc(data_lp.time{1}(1:10:end),[],data_lp.trial{1}(idx_meg,1:10:end)-repmat(mean(data_lp.trial{1}(idx_meg,1:10:end),2),1,size(data_lp.trial{1}(:,1:10:end),2)));
caxis([-1 1]*prctile(abs(dtmp(:)),99));
xlabel('time in s');
ylabel('channels');
title('filtered data, all channels whole exp.');

%% triggers
%  !! this part is totally idiosyncratic and will entirely depend on your
%  data / experiment. in general you can extract triggers using
%  event = ft_read_event(dsfile);
%  value = [event.value]';
%  sample = [event.sample]';
%
%  one common hickup: if the experiment uses button presses, these are
%  encoded as additional bits on the trigger channel, leading to the
%  addtition of 2^n to the trigger channel.
%  you can remove these and get clean triggers by calling dec2bin(value),
%  removing the additional lines, and calling bin2dec.
% 
%  what you need in the end:
%  - trialinfo, containing fields for each relevant variable
%  - trl, containing the sample of each stimulus / trial onset

load(bhvfile);

event = ft_read_event(dsfile);
idx_trig = find(strcmp({event.type},'UDIO001'));
value = [event(idx_trig).value];
sample = [event(idx_trig).sample];

value = value-256;

% identify trial numbers from preamble
idx_preamble_start = find(value==254);
idx_preamble_stop = find(value==255);

nwrong = [];
for trIDX=1:length(idx_preamble_start),
    trl_code = value(idx_preamble_start(trIDX)+[1:2]);
    trlnum(trIDX) = (trl_code(1)-1)*252+(trl_code(2)-1);

    if trIDX<length(idx_preamble_start),
        smps_act = sample(idx_preamble_stop(trIDX):idx_preamble_start(trIDX+1));
        trgs_act = value(idx_preamble_stop(trIDX):idx_preamble_start(trIDX+1));
    else
        smps_act = sample(idx_preamble_stop(trIDX):end);
        trgs_act = value(idx_preamble_stop(trIDX):end);
    end

    stimon_smps = smps_act(trgs_act==27);
    nstims = length(stimon_smps);
    n_bhv = sum(bhv.CodeNumbers{trlnum(trIDX)}==27);
    
    if n_bhv~=nstims,
        disp([num2str(trIDX) ' wrong n. bhv: ' num2str(n_bhv) ', trg: ' num2str(nstims)]);
        nwrong(trIDX) = 1;
    else
        nwrong(trIDX) = 0;
    end
  
    stimonset_trigs{trIDX} = stimon_smps;
    stim_coh{trIDX} = bhv.trialinfo.stims_coh{trlnum(trIDX)}(1:nstims);
    stim_ctr{trIDX} = bhv.trialinfo.stims_ctr{trlnum(trIDX)}(1:nstims);
    
end

keeptrials = find(~nwrong);
trl = cell2mat(stimonset_trigs(keeptrials));
trialinfo = [];
trialinfo.coh = cell2mat(stim_coh(keeptrials)');
trialinfo.ctr = cell2mat(stim_ctr(keeptrials)');

trl = round((trl./(hdr.Fs))*data_lp.fsample);

%% look at photodiode to get delay between trigger and stimulus
%  ideally, you record a photodiode signal to measure this delay yourself.
%  otherwise, the delay should usually be constant (if nothing goes wrong 
%  and the system is unchanged),so you could use the delay computed here 
%  (9ms) on your own data.

idx_uadc = find(cellfun(@(x) strcmp(x(1:3),'UAD'),data_lp.label));
idx_pd = idx_uadc(1);
dt_pd = data_hp.trial{1}(idx_pd,:);

dt_pd_aligned = [];
for trIDX=1:length(trl),
    dt_pd_aligned(trIDX,:) = dt_pd(trl(trIDX)+[-100:100]);
end
dt_pd_avg = mean(dt_pd_aligned);
h = hann(10);
dt_pd_avg = ms_convn(dt_pd_avg,h');
diff_pd = diff(dt_pd_avg);

figure;
plot([-100:100]./data_hp.fsample,[0 abs(diff(dt_pd_avg))]);
max_pre = prctile(abs(diff_pd(1:100)),90);
idx_first = find(abs(diff_pd)>(4*max_pre),1,'first');
t_delay = (idx_first-101)./data_hp.fsample;
hold on;
plot([t_delay t_delay],[-0.1 0.1],'k--');
xlabel('time from trigger');
title('derivative of photodiode response')

% correct onsets with the determined delay
trl = trl + round(t_delay*data_hp.fsample);


%% decoding using fs_cvmanova

% we downsample & filter the data further to make it faster and gain SNR

% lowpass filter
cfg_flt = [];
cfg_flt.lpfilter = 'yes';
cfg_flt.lpfreq = [30];
cfg_flt.lpinstabilityfix = 'split';
data_dec = ft_preprocessing(cfg_flt,data_lp);

% resample
cfg_rs = [];
cfg_rs.resamplefs = 100;
cfg_rs.detrend = 'no';
data_dec = ft_resampledata(cfg_rs,data_dec);

trl = round((trl./(data_lp.fsample))*data_dec.fsample);

prestim = round(data_dec.fsample*0.25);
poststim = round(data_dec.fsample*1);
trl = cat(2,trl'-prestim,trl'+poststim,-prestim*ones(size(trl')));

cfg_trial = [];
cfg_trial.trl = trl;
data_dec = ft_redefinetrial(cfg_trial,data_dec);

cfg_sel = [];
cfg_sel.channel = data_dec.label(idx_meg);
data_dec = ft_selectdata(cfg_sel,data_dec);

dat = reshape(cell2mat(data_dec.trial),length(data_dec.label),size(data_dec.trial{1},2),[]);

% check for noisy trials and channels. if the data is not too noisy it's
% usually not necessary to exclude anything (for decoding analyses useing
% low frequency broadband signals), but this would be a simple way
% to identify noisy data. the threshold used here is just a suggestion -
% depends on the data and noise level.

dat_var = squeeze(var(dat,[],2));
var_trials = mean(dat_var,1);
var_channels = mean(dat_var,2);

figure;
subplot(1,3,1);
imagesc(dat_var);
caxis([0 prctile(dat_var(:),95)])
xlabel('trials');
ylabel('channels');
title('variance');
subplot(1,3,2);
plot(zscore(var_trials));
hold on;
thresh_trials = prctile(zscore(var_trials),75)*10;
plot([1 length(var_trials)],[thresh_trials thresh_trials],'r');
xlabel('trials');
ylabel('z(var)');
title('variance over trials')
subplot(1,3,3);
plot(zscore(var_channels));
hold on;
thresh_chans = prctile(zscore(var_channels),75)*10;
plot([1 length(var_channels)],[thresh_chans thresh_chans],'r');
xlabel('channels');
ylabel('z(var)');
title('variance over channels')

% trials / channels to keep
idx_trl = find(zscore(var_trials)<thresh_trials);
idx_chn = find(zscore(var_channels)<thresh_chans);

% decode
cfg = [];
cfg.trialinfo = trialinfo;
cfg.trialinfo.coh = (cfg.trialinfo.coh(idx_trl)>2);
cfg.trialinfo.ctr = (cfg.trialinfo.ctr(idx_trl)>2);
cfg.varsDec = {'coh';'ctr'};
cfg.xvar = 1; % enable cross-variable decoding between coherence and contrast
cfg.nfold = 2;
cfg.tCov = 1:prestim;
cfg.avgCov = 1;

[D params] = fs_cvmanova(cfg,dat(idx_chn,:,idx_trl));

figure;
subplot(1,3,1);

% plot
cfg_plot = [];
ll=[];
cfg_plot.params = params;
cfg_plot.var = {'coh'};
cfg_plot.time = data_dec.time{1};
cfg_plot.color = [0.8 0.2 0.2];
ll(1)=fs_cvmanova_plot(cfg_plot,D);
hold on;
cfg_plot.var = {'ctr'};
cfg_plot.color = [0.2 0.2 0.8];
ll(2)=fs_cvmanova_plot(cfg_plot,D);
cfg_plot.var = {'ctr';'coh'};
cfg_plot.color = [0.2 0.8 0.2];
ll(3) = fs_cvmanova_plot(cfg_plot,D);

xlim([cfg_plot.time(1) cfg_plot.time(end)]);
plot([cfg_plot.time(1) cfg_plot.time(end)],[0 0],'k--');
legend(ll,{'coh';'ctr';'xvar'});
xlabel('time in s');
ylabel('neural information (D)');
title('within time');

% decode cross time
cfg = [];
cfg.trialinfo = trialinfo;
cfg.trialinfo.coh = cfg.trialinfo.coh(idx_trl);
cfg.trialinfo.ctr = cfg.trialinfo.ctr(idx_trl);
cfg.varsDec = {'coh';'ctr'};
cfg.xtime= 1; % enable cross time decoding
cfg.nfold = 2;
cfg.tCov = 1:prestim;
cfg.avgCov = 1;

[D params] = fs_cvmanova(cfg,dat(idx_chn,:,idx_trl));

subplot(1,3,2);
cfg_plot = [];
cfg_plot.params = params;
cfg_plot.var = {'coh'};
cfg_plot.time = data_dec.time{1};
cfg_plot.time2 = data_dec.time{1};
cfg_plot.xtime = 1;
fs_cvmanova_plot(cfg_plot,D);
caxis([-1.5 1.5]);
xlabel('time in s');
ylabel('time in s');
title('coherence cross time');
colorbar;
subplot(1,3,3);
cfg_plot.var = {'ctr'};
fs_cvmanova_plot(cfg_plot,D);
caxis([-1.5 1.5]);
xlabel('time in s');
ylabel('time in s');
title('contrast cross time');
colorbar;

%% source reconstruction

file_mri = [ftpath 'template/anatomy/single_subj_T1.nii'];
[vol, grad, grid, segmentedmri, leadfield, mri] = fs_coregister_meg_mri_leadfield(file_mri,dsfile);

% compute covariance matrix of data
C = double(reshape(dat,size(dat,1),[])*reshape(dat,size(dat,1),[])');

% compute filter
cfg_src = []; cfg_src.pow = 0; cfg_src.combine = 1;
source = fs_beamform(cfg_src,C,leadfield);

% apply filter to data (multiply them)
dat_long = reshape(dat,size(dat,1),[]);
dat_src_long = nan(length(source.pow),3,size(dat_long,2));
for sIDX=1:length(source.pow),
    dat_src_long(sIDX,:,:) =  squeeze(source.filt_raw(:,:,sIDX))*dat_long;
end

% data now has 3 directions per source (i.e. N_sources*3 channels)
label = {};
dat_src = [];
for dirIDX=1:size(dat_src_long,2),
    for srcIDX=1:size(dat_src_long,1),
        dat_src(end+1,:,:) = dat_src_long(srcIDX,dirIDX,:);
        label{end+1} = [num2str(srcIDX) '_' num2str(dirIDX)];
    end
end

dat_src = reshape(dat_src,size(dat_src,1),size(dat,2),size(dat,3));

% decode on the source level - whole brain
cfg = [];
cfg.trialinfo = trialinfo;
cfg.varsDec = {'coh';'ctr'};
cfg.nfold = 2;
cfg.tCov = 1:prestim;
cfg.avgCov = 1;
cfg.pca = 100; % we need to reduce dimensionality now, as we have too many channels!
cfg.pca_foldwise = 1;

[D_src params] = fs_cvmanova(cfg,dat_src);

figure;

% plot
cfg_plot = [];
ll=[];
cfg_plot.params = params;
cfg_plot.var = {'coh'};
cfg_plot.time = data_dec.time{1};
cfg_plot.color = [0.8 0.2 0.2];
ll(1)=fs_cvmanova_plot(cfg_plot,D_src);
hold on;
cfg_plot.var = {'ctr'};
cfg_plot.color = [0.2 0.2 0.8];
ll(2)=fs_cvmanova_plot(cfg_plot,D_src);

xlim([cfg_plot.time(1) cfg_plot.time(end)]);
plot([cfg_plot.time(1) cfg_plot.time(end)],[0 0],'k--');
legend(ll,{'coh';'ctr'});
xlabel('time in s');
ylabel('neural information (D)');
title('source level within time');


% define neighbours for searchlight analysis. for each source, we take the
% source itself and its direct neighbours in all 3 directions.
load('shell457.mat');
nbs = shell.neighbour;
for nbIDX=1:length(shell.neighbour), nbs{nbIDX}(end+1) = nbIDX; end;
for nbIDX=1:length(shell.neighbour), nbs{nbIDX} = cat(1,nbs{nbIDX},nbs{nbIDX}+457,nbs{nbIDX}+457*2); end;

% decode on the source level - searchlight
cfg = [];
cfg.trialinfo = trialinfo;
cfg.varsDec = {'coh';'ctr'};
cfg.nfold = 2;
cfg.tCov = 1:prestim;
cfg.avgCov = 1;
cfg.searchlight = 1;
cfg.neighbours = nbs;

[D params] = fs_cvmanova(cfg,dat_src);

% plot distribution of information on the source level
D = squeeze(mean(D,15));

figure;
subplot(1,2,1);
scatter3(shell.pos(:,1),shell.pos(:,2),shell.pos(:,3),80,squeeze(mean(D(:,1,:))),'filled');
title(cfg.varsDec{1});
axis equal
subplot(1,2,2);
scatter3(shell.pos(:,1),shell.pos(:,2),shell.pos(:,3),80,squeeze(mean(D(:,2,:))),'filled');
title(cfg.varsDec{2});
axis equal
