function [vol, grad, grid, segmentedmri, leadfield, mri] = fs_coregister_meg_mri_leadfield(file_mri,file_meg,varargin)
% sets fiducials in MRI scan, segments MRI and creates single shell
% headmodel. uses sensor positions from MEG data to create leadfield.
% inputs:
% file_mri - filename of the mri scan
% file_meg - filename of the raw meg data
% varargin{1} - optional file containing a source grid. default is a 457-point cortical shell.
% varargin{2} - optional gradiometer information to use. default is to read it out from the MEG header.
% varargin{3} - optional list of channels (cell array with channel names) to use.

if exist([file_mri '.mat']),
	load(file_mri);
else
	mri = ft_read_mri(file_mri);
end

% here you can select fiducials...
if ~isfield(mri,'coordsys') || ~strcmp(mri.coordsys,'ctf'),
	cfg = [];
	cfg.method = 'interactive';
	cfg.coordsys = 'ctf';
	mri = ft_volumerealign(cfg,mri);
end

% segment mri into skull, brain, ...
cfg           = [];
cfg.output    = 'brain';
segmentedmri  = ft_volumesegment(cfg, mri);
segmentedmri.anatomy = mri.anatomy;

% prepare 3D headmodel
cfg = [];
cfg.method='singleshell';
vol = ft_prepare_headmodel(cfg, segmentedmri);
vol = ft_convert_units(vol, 'cm');

hdr = ft_read_header(file_meg);

% optionally use list of channels
if length(varargin)>2 && ~isempty(varargin{3}),
    chan_use = varargin{3};
    ismeg = ismember(hdr.label,chan_use);
else
	% use all meg channels
    ismeg = ismember('meggrad',hdr.chantype);
end

% you can either submit a source grid as input to this function, or it uses
% by default the 457 shell
if length(varargin)>0 && ~isempty(varargin{1}),
	gridcfg = varargin{1};
else
	gridcfg = [];
	gridcfg.path = './shell457.mat';
	gridcfg.name = 'shell457';
	gridcfg.type = 'mat';
end

if strcmp(gridcfg.type,'mat'),
	load(gridcfg.path);
else
	shell = ft_read_headshape(gridcfg.path);
	shell.ori = normals(shell.pos,shell.tri,'vertex');
end

cfg                = [];
cfg.grid.pos       = shell.pos
cfg.grid.inside    = [1:size(cfg.grid.pos,1)];    % all of your ROIs are assumed to be inside the brain
cfg.grid.outside   = [];
if isfield(shell,'outside'),
	cfg.grid.inside = shell.inside;
	cfg.grid.outside = shell.outside;
end
cfg.grid.unit      = 'mm';
cfg.grid.nonlinear = 'yes';
cfg.grid.warpmni   = 'yes';
cfg.mri            = mri;

if ismeg,
    if length(varargin)>1 && ~isempty(varargin{2}),
        grad_use = varargin{2};
    else
        grad_use = hdr.grad;
    end
    cfg.grad = grad_use;
else
    % GET EEG ELEC POSITIONS!
end
cfg.template = cfg.grid;
grid = ft_prepare_sourcemodel(cfg); % template now in mri units!

if ~isfield(shell,'tri'),
	grid.ori = ones(size(grid.pos,1),3);
	grid.dist = pdist(grid.pos);
else
	grid.ori = normals(grid.pos,shell.tri,'vertex');
	grid.tri = shell.tri;
	
	if isfield(shell,'outside'),
		
		tmp = [];
		tmp.pnt = shell.ori(shell.outside,:);
		tmp = ft_transform_headshape(grid.params.Affine,tmp);
		tmp.pnt = tmp.pnt./repmat(sqrt(sum(tmp.pnt.^2,2)),1,3);
		
		grid.ori(shell.outside,:) = shell.ori(shell.outside,[2 1 3]);
		
    end
	
	if isfield(shell,'inside'),
		midlineidx = find((abs(shell.pos(shell.inside,1))<2.5).*((shell.pos(shell.inside,2)<30).*(shell.pos(shell.inside,2)>-45)));
		grid.dist = fs_mesh_dist(grid.pos(shell.inside,:),grid.tri,midlineidx);
	else
		midlineidx = find((abs(shell.pos(:,1))<2.5).*((shell.pos(:,2)<30).*(shell.pos(:,2)>-45)));
		grid.dist = fs_mesh_dist(grid.pos,grid.tri,midlineidx);
	end
end

grid2 = ft_convert_units(grid,'cm');
grid2.params = grid.params;
grid2.cfg = grid.cfg;

if ~isfield(grid2,'outside'),
    grid2.outside = logical(zeros(size(grid2.inside)));
end

figure;
ft_plot_headshape(grid2)
%ft_plot_mesh(grid.pos(grid.inside,:)./10);
hold on
ft_plot_headmodel(vol, 'edgecolor', 'none','facecolor',[0.2 0.8 0.2]); alpha 0.5;
ft_plot_sens(grad_use);

cfg                 = [];
cfg.grad            = grad_use;
cfg.vol             = vol;
cfg.channel         = {'MEG'};
cfg.grid            = grid2;

insidetmp = cfg.grid.inside;
if islogical(cfg.grid.outside),
    outsidetmp = find(cfg.grid.outside);
else
    outsidetmp = cfg.grid.outside;
end

cfg.grid.inside = 1:size(cfg.grid.pos,1);
cfg.grid.outside = [];

[leadfield] = ft_prepare_leadfield(cfg);

leadfield.inside(outsidetmp) = 0;
if ~isfield(leadfield,'outside'),
leadfield.outside = logical(zeros(size(leadfield.inside)));
end
grad = grad_use;

end

        
       
   