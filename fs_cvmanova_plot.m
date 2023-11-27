function handles=fs_cvmanova_plot(cfg,D),

%% inputs:
% cfg.var:    1x1 cell array which variable to plot:
%             cfg.var = {'stim'};
%             for cross variable decoding: 2x1 cell array:
%             cfg.var = {'stim';'resp'};
% cfg.split:  cell array specifying which variables to split by, and which
%             level of the split to plot. e.g. if you want to plot
%             information at stim=2 and resp=3:
%             cfg.split = {{'stim';2};{'resp';3}};
% cfg.split2: for cross-split decoding specify the test split (as in cfg.split)
%             cfg.split is then the training split
% cfg.xtime:  whether to plot cross time information (as an imagesc plot)
% cfg.time:   time vector to plot with accurate time axis
% cfg.time2:  time vector to plot with accurate time axis (2nd time
%             dimension for asymmetric matrices)
% cfg.color:  3x1 vector containing rgb values
% cfg.
% 
% D:          [time  x split1   x ... x split5   x vars
%             time2 x split2_1 x ... x split2_5 x vars2 x folds x subjects] (16D array)
%             "pattern distinctness", i.e. roughly variability explained by a
%             certain contrast in relation to noise variability. for all splitX,
%             the first index contains the overall effect, the subsequent ones
%             contain effects when considering trials with the n-1st value of the
%             split variable.
%             to assemble an "all-subjects"-D from all the single-subject
%             Ds: D_new = cat(16,D1,D2,D3,...,Dn);
%             examples:
%             A) 10 time points, 2 variables, 0 splits, 5 folds, 4 subjects, no cross-decoding.
%                D = [10 x 1 x 1 x 1 x 1 x 1 x 2 x 1 x 1 x 1 x 1 x 1 x 1 x 1 x 5 x 4]
%             B) same, but cross-variable decoding
%                D = [10 x 1 x 1 x 1 x 1 x 1 x 2 x 1 x 1 x 1 x 1 x 1 x 1 x 2 x 5 x 4]
%             C) same, but also cross-time decoding
%                D = [10 x 1 x 1 x 1 x 1 x 1 x 2 x 10 x 1 x 1 x 1 x 1 x 1 x 2 x 5 x 4]
%             D) same, but splitting across one variable and doing
%                cross-split decoding
%                D = [10 x 3 x 1 x 1 x 1 x 1 x 2 x 10 x 3 x 1 x 1 x 1 x 1 x 2 x 5 x 4]
%             

if ~isfield(cfg,'var'),  error('specify variable to plot!'); end
if ~isfield(cfg,'xtime'),  cfg.xtime = 0; end
if ~isfield(cfg,'split'),  cfg.split = {}; end
if ~isfield(cfg,'split2'),  cfg.split2 = {}; end % for cross-split: train on split, test on split2

if ~isfield(cfg,'time'), cfg.time = 1:size(D,1); end
if ~isfield(cfg,'time2'), cfg.time2 = 1:size(D,8); end
if ~isfield(cfg,'color'), cfg.color = [0 0 0]; end

if ~isfield(cfg,'transparent'), cfg.transparent = 0; end
if ~isfield(cfg,'saturation'), cfg.saturation = 0.15; end

if ~isfield(cfg,'stats'), cfg.stats = []; end


%if ~isfield(cfg,'xvar'),   cfg.xvar = 0; end % redundant!!
if length(cfg.var)==2, cfg.xvar=1; else cfg.xvar = 0; end
%if ~isfield(cfg,'xsplit'), cfg.xsplit = 0; end % redundant!!
if ~isempty(cfg.split2), cfg.xsplit = 1; else cfg.xsplit = 0; end

vars = cfg.params.cfg.varsDec;
varsSplit = cfg.params.cfg.varsSplit;


% which variable?
var_idx = find(strcmp(vars,cfg.var{1}));

if cfg.params.cfg.xvar==0, var_idx2 = 1;
elseif cfg.xvar, var_idx2 = find(strcmp(vars,cfg.var{2}));
else var_idx2 = var_idx; end

% which split
if isempty(cfg.split),
    split_idx  = [1 1 1 1 1];
    split_idx2 = [1 1 1 1 1];
else
    % split...
    split_idx  = [1 1 1 1 1];
    split_idx2 = [1 1 1 1 1];
    
    for i=1:length(cfg.split),
        split_idx(find(strcmp(varsSplit,cfg.split{i}{1}))) = cfg.split{i}{2}+1;
    end
    
    if cfg.xsplit,
        for i=1:length(cfg.split2),
            split_idx2(find(strcmp(varsSplit,cfg.split2{i}{1}))) = cfg.split2{i}{2}+1;
        end
	else
		if cfg.params.cfg.xsplit,
			split_idx2 = split_idx;
		else
			
		end
    end
end

Dplot = D(:,split_idx(1),split_idx(2),split_idx(3),split_idx(4),split_idx(5),var_idx,:,split_idx2(1),split_idx2(2),split_idx2(3),split_idx2(4),split_idx2(5),var_idx2,:,:,:);
Dplot = mean(mean(Dplot,15),17);


if cfg.xtime,
    handles = imagesc(cfg.time2,cfg.time,squeeze(mean(Dplot,16)));
else
    
    if cfg.params.cfg.xtime,
        
        idxdiag = eye(size(Dplot,1));
        
        idxdiag = permute(idxdiag,[1 3 4 5 6 7 8 2 9 10 11 12 13 14 15 16]);
        idxdiag = repmat(idxdiag,[1 size(Dplot,2) size(Dplot,3) size(Dplot,4) size(Dplot,5) size(Dplot,6) size(Dplot,7) 1 size(Dplot,9) size(Dplot,10) size(Dplot,11) size(Dplot,12) size(Dplot,13) size(Dplot,14) size(Dplot,15) size(Dplot,16)]);
        
        Dplot = reshape(Dplot(find(idxdiag)),[size(Dplot,1) size(Dplot,2) size(Dplot,3) size(Dplot,4) size(Dplot,5) size(Dplot,6) size(Dplot,7) 1 size(Dplot,9) size(Dplot,10) size(Dplot,11) size(Dplot,12) size(Dplot,13) size(Dplot,14) size(Dplot,15) size(Dplot,16)]);
    end
    
    if size(Dplot,16)==1,
        handles = plot(cfg.time,squeeze(mean(Dplot,16)),'Color',cfg.color,'LineWidth',2);
    else
        handles = shadedErrorBar(cfg.time,squeeze(mean(Dplot,16)),squeeze(std(Dplot,[],16))./sqrt(sum(~isnan(Dplot(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,:)))),{'Color',cfg.color,'LineWidth',2},cfg.transparent,cfg.saturation);
    end
    
    
end




end

function varargout=shadedErrorBar(x,y,errBar,lineProps,transparent,saturation)
% function H=shadedErrorBar(x,y,errBar,lineProps,transparent)
%
% Purpose
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
% errBar - if a vector we draw symmetric errorbars. If it has a size
%          of [2,length(x)] then we draw asymmetric error bars with
%          row 1 being the upper bar and row 2 being the lower bar
%          (with respect to y). ** alternatively ** errBar can be a
%          cellArray of two function handles. The first defines which
%          statistic the line should be and the second defines the
%          error bar.
% lineProps - [optional,'-k' by default] defines the properties of
%             the data line. e.g.:
%             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image.
%
% Outputs
% H - a structure of handles to the generated plot objects.
%
%
% Examples
% y=randn(30,80); x=1:size(y,2);
% shadedErrorBar(x,mean(y,1),std(y),'g');
% shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});
% shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'});
%
% Overlay two transparent lines
% y=randn(30,80)*10; x=(1:size(y,2))-40;
% shadedErrorBar(x,y,{@mean,@std},'-r',1);
% hold on
% y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
% shadedErrorBar(x,y,{@mean,@std},'-b',1);
% hold off
%
%
% Rob Campbell - November 2009



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error checking
error(nargchk(3,6,nargin))


%Process y using function handles if needed to make the error bar
%dynamically
if iscell(errBar)
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(y);
    y=fun1(y);
else
    y=y(:)';
end

if isempty(x)
    x=1:length(y);
else
    x=x(:)';
end


%Make upper and lower error bars if only one was specified
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    s=size(errBar);
    f=find(s==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(x) ~= length(errBar)
    error('length(x) must equal length(errBar)')
end

%Set default options
defaultProps={'-k'};
if nargin<4, lineProps=defaultProps; end
if isempty(lineProps), lineProps=defaultProps; end
if ~iscell(lineProps), lineProps={lineProps}; end

if nargin<5, transparent=0; end

if nargin<6, saturation=0.15; end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot to get the parameters of the line
H.mainLine=plot(x,y,lineProps{:});


% Work out the color of the shaded region and associated lines
% Using alpha requires the render to be openGL and so you can't
% save a vector image. On the other hand, you need alpha if you're
% overlaying lines. There we have the option of choosing alpha or a
% de-saturated solid colour for the patch surface .

col=get(H.mainLine,'color');
edgeColor=col+(1-col)*0.55;
patchSaturation=saturation; %How de-saturated or transparent to make patch
if transparent
    faceAlpha=patchSaturation;
    patchColor=col;
    set(gcf,'renderer','openGL')
else
    faceAlpha=1;
    patchColor=col+(1-col)*(1-patchSaturation);
    set(gcf,'renderer','painters')
end


%Calculate the error bars
uE=y+errBar(1,:);
lE=y-errBar(2,:);


%Add the patch error bar
holdStatus=ishold;
if ~holdStatus, hold on,  end


%Make the patch
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];

%remove nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];


H.patch=patch(xP,yP,1,'facecolor',patchColor,...
    'edgecolor','none',...
    'facealpha',faceAlpha);


%Make pretty edges around the patch.
%H.edge(1)=plot(x,lE,'-','color',edgeColor);
%H.edge(2)=plot(x,uE,'-','color',edgeColor);

%Now replace the line (this avoids having to bugger about with z coordinates)
delete(H.mainLine)
H.mainLine=plot(x,y,lineProps{:});


if ~holdStatus, hold off, end


if nargout==1
    varargout{1}=H;
end

end

