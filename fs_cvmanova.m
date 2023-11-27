function [D params] = fs_cvmanova(cfg,data),

if ~isfield(cfg,'xtime'), cfg.xtime = 0; end
if ~isfield(cfg,'xvar'), cfg.xvar = 0; end
if ~isfield(cfg,'xsplit'), cfg.xsplit = 0; end
if ~isfield(cfg,'nfold'), cfg.nfold = 10; end
if ~isfield(cfg,'tCov'), cfg.tCov = 1:size(data,2); end
if ~isfield(cfg,'lambda'), cfg.lambda = 0; end
if ~isfield(cfg,'varsSplit'), cfg.varsSplit = {}; end
if ~isfield(cfg,'varsDec'), cfg.varsDec = fieldnames(cfg.trialinfo); end
if ~isfield(cfg,'groundTruth'), cfg.groundTruth = 0; end
if ~isfield(cfg,'single'), cfg.single = 0; end
if ~isfield(cfg,'searchlight'), cfg.searchlight = 0; end
if ~isfield(cfg,'neighbours'), cfg.neighbours = {1:size(data,1)}; end
if ~isfield(cfg,'singletrial'), cfg.singletrial = 0; end
if ~isfield(cfg,'avgCov'), cfg.avgCov = 0; end
if ~isfield(cfg,'disableOutput'), cfg.disableOutput = 0; end
if ~isfield(cfg,'cv'), cfg.cv = []; end
if ~isfield(cfg,'saveE'), cfg.saveE = 0; end
if ~isfield(cfg,'paramsE'), cfg.paramsE = []; end
if ~isfield(cfg,'pca'), cfg.pca = 0; end
if ~isfield(cfg,'pca_average'), cfg.pca_average = 1; end
if ~isfield(cfg,'pca_foldwise'), cfg.pca_foldwise = 0; end

if isfield(cfg,'rndseed'),
    rng(cfg.rndseed);
end


if cfg.single,
    data = single(data);
end

trialinfo = [];
for vIDX=1:length(cfg.varsDec),
    if cfg.single,
        eval(['trialinfo.' cfg.varsDec{vIDX} '=single(cfg.trialinfo.' cfg.varsDec{vIDX} ');']);
    else
        eval(['trialinfo.' cfg.varsDec{vIDX} '=cfg.trialinfo.' cfg.varsDec{vIDX} ';']);
    end
end

[X sid cs sc sm] = make_dm([],trialinfo);
cs = cellfun(@(x) x',cs,'UniformOutput',0);

if cfg.pca,
    
    fE = floor((size(data,3)/cfg.nfold)*(cfg.nfold-1))-rank(X);
    maxDim = (fE-1)-1;
    
    if cfg.pca>1,
        maxDim = min(cfg.pca,size(data,1));
    else
        maxDim = min([maxDim size(data,1)]);
    end
    
    if ~cfg.pca_foldwise,
        if cfg.pca_average,
            dtm = [];
            for cndIDX=1:max(sid),
                
                dtm(1:size(data,1),:,cndIDX) = mean(data(:,:,sid==cndIDX),3);
            end
            dtm = reshape(dtm,size(dtm,1),[]);
        else
            dtm = reshape(data,size(data,1),[]);
            dtm = dtm(:,randperm(size(dtm,2),min(size(dtm,2),1000000)));
        end
        [coeff score latent] = pca(dtm');
        dt = reshape(data,size(data,1),[]);
        ddim = size(data);
        if ddim(1)>size(coeff,2),
            ddim(1) = size(coeff,2);
        end
        dt = dt'*coeff;
        dt = reshape(dt',ddim);
        
        dt = dt(1:maxDim,:,:);
        
        data = dt;
        
        clear dt
        
        prc_var = sum(latent(1:maxDim))./sum(latent);
    end
else
    prc_var = 1;
end

conds_split = {};
for i=1:5,
	conds_split{i}{1} = ones(1,size(cs{1},2));
end

for vsIDX=1:length(cfg.varsSplit),
	idx_splitvar = find(strcmp(cfg.varsDec,cfg.varsSplit{vsIDX}));
	splitvals = unique(cs{idx_splitvar});
	conds_split{vsIDX}{1} = ones(1,size(cs{1},2));
	for svIDX = 1:length(splitvals),
		conds_split{vsIDX}{svIDX+1} = cs{idx_splitvar}==splitvals(svIDX);
	end
end




if isempty(cfg.cv),
    if ~cfg.groundTruth,
		if cfg.nfold==1,
			cv = [];
			cv.NumTestSets = 1;
			cv.NumObservations = size(data,3);
			cv.training{1} = logical(ones(1,size(data,3)));
			cv.test{1} = logical(ones(1,size(data,3)));
        else
            
            cv = cvpartition(sid,'Kfold',cfg.nfold);
		end
        
    else
        cv = [];
        cv.NumTestSets = 1;
        cv.training(1,:) = logical(ones(1,size(data,3)));
        cv.test(1,:) = logical(ones(1,size(data,3)));
    end
else
    cv = cfg.cv;
end


if cfg.searchlight,
	nSL = length(cfg.neighbours);
else
	nSL = 1;
	cfg.neighbours = {1:size(data,1)};
end


nT = size(data,2);
if cfg.xtime, nT2 = nT; else nT2 = 1; end
nV = length(cs);
if cfg.xvar, nV2 = nV; else nV2 = 1; end
nS11 = length(conds_split{1});
nS21 = length(conds_split{2});
nS31 = length(conds_split{3});
nS41 = length(conds_split{4});
nS51 = length(conds_split{5});
if cfg.xsplit,
	nS12 = nS11; nS22 = nS21; nS32 = nS31; nS42 = nS41; nS52 = nS51;
else
	nS12 = 1; nS22 = 1; nS32 = 1; nS42 = 1; nS52 = 1;
end

if cfg.singletrial,
	try
		D = nan(nT,nS11,nS21,nS31,nS41,nS51,nV,nT2,nS12,nS22,nS32,nS42,nS52,nV2,cv.NumObservations,nSL);
	catch
		D = nan(nT,nS11,nS21,nS31,nS41,nS51,nV,nT2,nS12,nS22,nS32,nS42,nS52,nV2,cv.N,nSL);
	end
else
    D = nan(nT,nS11,nS21,nS31,nS41,nS51,nV,nT2,nS12,nS22,nS32,nS42,nS52,nV2,cv.NumTestSets,nSL);
end
if cfg.single,
    D = single(D);
end


for slIDX=1:nSL,
	
data_sl = data(cfg.neighbours{slIDX},:,:);
if cfg.disableOutput<3,
    disp(['sl' num2str(slIDX) ' of ' num2str(nSL)]);
end
nChan = size(data_sl,1);
params_all = {};
params_mean = {};
params_test = {};

for kIDX = 1:cv.NumTestSets,
                    
    if strcmp(class(cv),'cvpartition'),
        idxtrain = cv.training(kIDX);
        idxtest = cv.test(kIDX);
    else
        idxtrain = cv.training{kIDX};
        idxtest = cv.test{kIDX};
    end
    
    data_train = data(:,:,idxtrain);
   
    if cfg.pca && cfg.pca_foldwise,
        [data_train coeff{kIDX} latent{kIDX} prc_var(kIDX)] = do_pca(cfg.pca_average,data_train,sid(idxtrain),maxDim);
    end
    
end

if isempty(cfg.paramsE),
    
    if ~isstr(cfg.tCov),
        t0 = tic;
        if cfg.disableOutput<2,
            disp('computing covariance. ');
        end
        if ~cfg.groundTruth,
            params_all = {};
			if cfg.avgCov,
				params_mean = {};
				for kIDX = 1:cv.NumTestSets,
                    
                    if strcmp(class(cv),'cvpartition'),
                        idxtrain = cv.training(kIDX);
                        idxtest = cv.test(kIDX);
                    else
                        idxtrain = cv.training{kIDX};
                        idxtest = cv.test{kIDX};
                    end
                    
                    data_act = squeeze(mean(data_sl(:,cfg.tCov,:),2));
                    if size(data_act,2)==1, data_act = data_act'; end
                    
                    data_train = data_act(:,idxtrain)';
                    data_test  = data_act(:,idxtest)';
                    X_train = X(idxtrain,:);
                    X_test = X(idxtest,:);
                    
                    if cfg.pca && cfg.pca_foldwise,
                        [data_train] = do_pca(cfg.pca_average,data_train',sid(idxtrain),maxDim,coeff{kIDX},latent{kIDX});
                        [data_test] = do_pca(cfg.pca_average,data_test',sid(idxtest),maxDim,coeff{kIDX},latent{kIDX});
                        data_train = data_train';
                        data_test = data_test';
                    end
                    
					[params_mean{kIDX}] = fs_manova_train(data_train,X_train,cs,'lambda',cfg.lambda);
				end
			else
				for tIDX=1:length(cfg.tCov),
					for kIDX = 1:cv.NumTestSets,
                        
                        
                        if strcmp(class(cv),'cvpartition'),
                            idxtrain = cv.training(kIDX);
                            idxtest = cv.test(kIDX);
                        else
                            idxtrain = cv.training{kIDX};
                            idxtest = cv.test{kIDX};
                        end
                        
						data_act = reshape(squeeze(data_sl(:,cfg.tCov(tIDX),:)),[],size(data_sl,3)); 
						if size(data_act,2)==1, data_act = data_act'; end
						
                        data_train = data_act(:,idxtrain)';
                        
                        data_test  = data_act(:,idxtest)';
                        X_train = X(idxtrain,:);
                        X_test = X(idxtest,:);
                        
                        if cfg.pca && cfg.pca_foldwise,
                            [data_train] = do_pca([],data_train',[],maxDim,coeff{kIDX},latent{kIDX});
                            [data_test] = do_pca([],data_test',[],maxDim,coeff{kIDX},latent{kIDX});
                            data_train = data_train';
                            data_test = data_test';
                        end
                        
						[params_all{tIDX,kIDX}] = fs_manova_train(data_train,X_train,cs,'lambda',cfg.lambda);
					end
				end
				
                params_mean = {};
                for tIDX=1:length(cfg.tCov),
                    for kIDX = 1:cv.NumTestSets,
                        if tIDX==1,
                            params_mean{kIDX}.iE =  params_all{tIDX,kIDX}.iE./length(cfg.tCov);
                            params_mean{kIDX}.E =  params_all{tIDX,kIDX}.E./length(cfg.tCov);
                        else
                            params_mean{kIDX}.iE =  params_mean{kIDX}.iE+params_all{tIDX,kIDX}.iE./length(cfg.tCov);
                            params_mean{kIDX}.E =  params_mean{kIDX}.E+params_all{tIDX,kIDX}.E./length(cfg.tCov);
                        end
                        params_mean{kIDX}.fE = params_all{tIDX,kIDX}.fE;
                    end
                end
                
            end
        else
            params_mean{1}.iE = cfg.iE;
            params_mean{1}.fE = 0;
        end
        t1 = toc(t0);
        if cfg.disableOutput<2,
            fprintf('\b%s',[' time elapsed: ' num2str(t1)]);
        end
    end
else
    for kIDX = 1:cv.NumTestSets,
        params_mean{kIDX}.iE = cfg.paramsE.iE{kIDX};
        params_mean{kIDX}.fE = cfg.paramsE.fE{kIDX};
    end
end

nTStr = sprintf('%04d',nT);
t0 = tic;
timStr = sprintf('%06d',0);

if cfg.disableOutput<2,
    fprintf('\n%s',['training 0000 of ' nTStr '. time elapsed: ' timStr 's.']);
end
for tIDX1 = 1:nT,
	t1 = round(toc(t0));
	timStr = sprintf('%06d',t1);
	TStr = sprintf('%04d',tIDX1);

    if cfg.disableOutput<1,
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s',[TStr ' of ' nTStr '. time elapsed: ' timStr 's.']);
    end
    
	for s1IDX1 = 1:length(conds_split{1}),
	for s2IDX1 = 1:length(conds_split{2}),
	for s3IDX1 = 1:length(conds_split{3}),
	for s4IDX1 = 1:length(conds_split{4}),
	for s5IDX1 = 1:length(conds_split{5}),

        conds_use = conds_split{1}{s1IDX1}.*conds_split{2}{s2IDX1}.*conds_split{3}{s3IDX1}.*conds_split{4}{s4IDX1}.*conds_split{5}{s5IDX1};
			
		cs_use = {};
		for csIDX=1:length(cs),
			cs_use{csIDX} = cs{csIDX};
            cs_use{csIDX}(:,find(~conds_use)) = 0;
		end
		
		for kIDX = 1:cv.NumTestSets,
			
			if ~cfg.groundTruth,
				if strcmp(class(cv),'cvpartition')
					data_train = reshape(data_sl(:,tIDX1,cv.training(kIDX)),nChan,sum(cv.training(kIDX)))';
					X_train = X(cv.training(kIDX),:);
				else
					data_train = reshape(data_sl(:,tIDX1,cv.training{kIDX}),nChan,sum(cv.training{kIDX}))';
					X_train = X(cv.training{kIDX},:);
                end
                
                if cfg.pca && cfg.pca_foldwise,
                    [data_train] = do_pca([],data_train',[],maxDim,coeff{kIDX},latent{kIDX});
                    data_train = data_train';
                    
                end
                
				[params_train{tIDX1,kIDX,s1IDX1,s2IDX1,s3IDX1,s4IDX1,s5IDX1}] = fs_manova_train(data_train,X_train,cs_use,'lambda',cfg.lambda);
			else
				data_train = reshape(data_sl(:,tIDX1,cv.training(kIDX)),nChan,sum(cv.training(kIDX)))'; 
				X_train = X(cv.training(kIDX,:),:);
                
                if cfg.pca && cfg.pca_foldwise,
                    [data_train] = do_pca([],data_train',[],maxDim,coeff{kIDX},latent{kIDX});
                    data_train = data_train';
                           
                end
              
				[params_train{tIDX1,kIDX,s1IDX1,s2IDX1,s3IDX1,s4IDX1,s5IDX1}] = fs_manova_train(data_train,X_train,cs_use,'lambda',cfg.lambda,'beta',squeeze(cfg.betaTrue(find(conds_use),:,tIDX1)));
            end
            
            
            
		end
	end
	end
	end
	end
	end
	
end

t0 = tic;
timStr = sprintf('%06d',0);
if cfg.disableOutput<2,
    fprintf('\n%s',['testing  0000 of ' nTStr '. time elapsed: ' timStr 's.']);
end
for tIDX1 = 1:nT,
	
	t1 = round(toc(t0));
	timStr = sprintf('%06d',t1);
	TStr = sprintf('%04d',tIDX1);
	
    if cfg.disableOutput<1,
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%s',[TStr ' of ' nTStr '. time elapsed: ' timStr 's.']);
    end
    
	for s1IDX1 = 1:length(conds_split{1}),
	for s2IDX1 = 1:length(conds_split{2}),
	for s3IDX1 = 1:length(conds_split{3}),
	for s4IDX1 = 1:length(conds_split{4}),
	for s5IDX1 = 1:length(conds_split{5}),
		
		for tIDX2 = 1:nT2,
			
			for s1IDX2 = 1:nS12,
			for s2IDX2 = 1:nS22,
			for s3IDX2 = 1:nS32,
			for s4IDX2 = 1:nS42,
			for s5IDX2 = 1:nS52,
				
				sAIDX1 = [s1IDX1 s2IDX1 s3IDX1 s4IDX1 s5IDX1];
				sAIDX2 = [s1IDX2 s2IDX2 s3IDX2 s4IDX2 s5IDX2];
				
				if cfg.xsplit,
					badsplit = sum((sAIDX2==1).*(sAIDX1~=1))>0;
					subsplit = sum((sAIDX1==1).*(sAIDX2~=1))>0;
				else
					badsplit = 0;
					subsplit = 0;
				end
			
				if ~badsplit,
					
					if subsplit,
						a=1;
					end
					
					conds_use_train = conds_split{1}{s1IDX1}.*conds_split{2}{s2IDX1}.*conds_split{3}{s3IDX1}.*conds_split{4}{s4IDX1}.*conds_split{5}{s5IDX1};
					
					if cfg.xsplit,
						conds_use = conds_split{1}{s1IDX2}.*conds_split{2}{s2IDX2}.*conds_split{3}{s3IDX2}.*conds_split{4}{s4IDX2}.*conds_split{5}{s5IDX2};
						s1IDX2act = s1IDX2;
						s2IDX2act = s2IDX2;
						s3IDX2act = s3IDX2;
						s4IDX2act = s4IDX2;
						s5IDX2act = s5IDX2;
					else
						conds_use = conds_use_train;
						s1IDX2act = s1IDX1;
						s2IDX2act = s2IDX1;
						s3IDX2act = s3IDX1;
						s4IDX2act = s4IDX1;
						s5IDX2act = s5IDX1;
					end
					
					if cfg.xtime,
						tIDX = tIDX2;
					else
						tIDX = tIDX1;
					end
					
					cs_use = {};
					for csIDX=1:length(cs),
						cs_use{csIDX} = cs{csIDX};
                        cs_use{csIDX}(:,find(~conds_use)) = 0;
					end
					
					cs_use_train = {};
					for csIDX=1:length(cs),
						cs_use_train{csIDX} = cs{csIDX};
                        cs_use_train{csIDX}(:,find(~conds_use_train)) = 0;
					end
					
					for kIDX = 1:cv.NumTestSets,
						
						if ~cfg.groundTruth,
							if strcmp(class(cv),'cvpartition'),
								data_test = reshape(data_sl(:,tIDX,cv.test(kIDX)),nChan,sum(cv.test(kIDX)))';
								X_test= X(cv.test(kIDX),:);
							else
								data_test = reshape(data_sl(:,tIDX,cv.test{kIDX}),nChan,sum(cv.test{kIDX}))';
								X_test= X(cv.test{kIDX},:);
							end
						else
							 data_test = reshape(data_sl(:,tIDX,cv.test(kIDX)),nChan,sum(cv.test(kIDX)))';
							X_test= X(cv.test(kIDX,:),:);
                        end
                       
                        if cfg.pca && cfg.pca_foldwise,
                            [data_test] = do_pca([],data_test',[],maxDim,coeff{kIDX},latent{kIDX});
                          
                            data_test = data_test';
                        end
                      
                        if ~isempty(cfg.correct),
                            if strcmp(class(cv),'cvpartition'),
                                idxtest = cv.test(kIDX);
                            else
                                idxtest = cv.test{kIDX};
                            end
                            [data_test avg_correct] = correct_vars(data_test,cfg.correct(idxtest,:),[],cfg.correct_avg);
                        end
                        
						params_act = params_train{tIDX1,kIDX,s1IDX1,s2IDX1,s3IDX1,s4IDX1,s5IDX1};
						params_act.iE = params_mean{kIDX}.iE;
						params_act.fE = params_mean{kIDX}.fE;
						
						params_act_pre = params_act;
						
						for v1IDX = 1:nV2,
							for v2IDX = 1:nV2,
								
								params_act = params_act_pre;
								
								if nV2==1, 
									X_test_act = X_test;
                                    for v1IDXB = 1:nV,
                                        cs_var1 = cs_use_train{v1IDXB};
                                        cs_var2 = cs_use{v1IDXB};
                                        csfct = sum(cs_var2~=0)/sum(cs_var1~=0);
                                        params_act.CCs{v1IDXB} = (pinv(cs_var1)*cs_var2)/csfct;
                                    end  
								else
									cs_var1 = cs_use_train{v1IDX};
									cs_var2 = cs_use{v2IDX};
                                    
                                	X_test_act = X_test;
                                    csfct = sum(cs_var2~=0)/sum(cs_var1~=0);
                                    params_act.CCs = {(pinv(cs_var1)*cs_var2)/csfct};
                                  	params_act.betaDelta = params_act.betaDelta(v1IDX);
									
                                end
                                
                                try
                                    params_act_test = params_test{tIDX,kIDX,s1IDX2act,s2IDX2act,s3IDX2act,s4IDX2act,s5IDX2act,v1IDX,v2IDX};
                                    if ~isempty(params_act_test),
                                        params_act.beta_test = params_act_test.beta_test;
                                        params_act.XX = params_act_test.XX;
                                        params_act.n_used = params_act_test.n_used;
                                   end
                                catch
                                    
                                end
                                
                                try
                                    if size(X_test_act,1)==size(data_test,1),
                                        if ~cfg.groundTruth,
                                            if ~cfg.singletrial,
                                                [Dact, ~,  params_test{tIDX,kIDX,s1IDX2act,s2IDX2act,s3IDX2act,s4IDX2act,s5IDX2act,v1IDX,v2IDX}]          = fs_manova_test(params_act,data_test,X_test_act);
                                            else
                                                params_act.singleTrial = 1;
                                                [~, ~,  params_test{tIDX,kIDX,s1IDX2act,s2IDX2act,s3IDX2act,s4IDX2act,s5IDX2act,v1IDX,v2IDX},Dact] = fs_manova_test(params_act,data_test,X_test_act);
                                            end
                                        else
                                            beta_test = squeeze(cfg.betaTrue(find(conds_use),:,tIDX));
                                            params_act.cfact = 1/sum(X_test_act(:));
                                            if nV2==1,
                                                params_act.beta_test = beta_test;
                                            else
                                                params_act.beta_test(v11,:) = beta_test(v21,:);
                                                params_act.beta_test(v1n1,:) = beta_test(v2n1,:);
                                            end
                                            [Dact, ~,  params_test{tIDX,kIDX,s1IDX2act,s2IDX2act,s3IDX2act,s4IDX2act,s5IDX2act,v1IDX,v2IDX}] = fs_manova_test(params_act,data_test,X_test_act);
                                        end
                                    else
                                        Dact = NaN;
                                    end
                                catch ME
									Dact = NaN;
                                end
								
                                if cfg.singletrial,
									if strcmp(class(cv),'cvpartition'),
										foldIDS = find(cv.test(kIDX));
									else
										 foldIDS = find(cv.test{kIDX});
									end
                                else
                                    foldIDS = kIDX;
                                end
                                
								if nV2==1,
                                    D(tIDX1,s1IDX1,s2IDX1,s3IDX1,s4IDX1,s5IDX1,:,tIDX2,s1IDX2,s2IDX2,s3IDX2,s4IDX2,s5IDX2,1,foldIDS,slIDX) = Dact;
								else
									D(tIDX1,s1IDX1,s2IDX1,s3IDX1,s4IDX1,s5IDX1,v1IDX,tIDX2,s1IDX2,s2IDX2,s3IDX2,s4IDX2,s5IDX2,v2IDX,foldIDS,slIDX) = Dact;
								end
								
							end
						end
						
					end
					
				else
					D(tIDX1,s1IDX1,s2IDX1,s3IDX1,s4IDX1,s5IDX1,:,tIDX2,s1IDX2,s2IDX2,s3IDX2,s4IDX2,s5IDX2,:,:,slIDX) = NaN;
				end
				
			end
			end
			end
			end
			end
			
		end
		
	end
	end
	end
	end
	end
end

params = [];
params.cv = cv;
params.cfg = cfg;
params.X   = X;
params.sid = sid;
params.cs = cs;
params.sc = sc;
params.sm = sm;

params.var = prc_var;

if cfg.saveE,
    for kIDX=1:cv.NumTestSets,
        params.iE{kIDX} = params_mean{kIDX}.iE;
        params.fE{kIDX} = params_mean{kIDX}.fE;
    end
end

if cfg.disableOutput<2,
    fprintf('\n');
end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEPENDENCIES

function [dm stratid cs stratcombs stratmat] = make_dm(cfg,trialinfo),

if ~isfield(cfg,'vars'), cfg.vars = fieldnames(trialinfo); end

for vIDX=1:length(cfg.vars),
	labels(:,vIDX) = eval(['trialinfo.' cfg.vars{vIDX} ';']);
end

ll = unique(labels(:));
ll = ll(~isnan(ll));
nanval = max(ll)+1;

labels(isnan(labels(:))) = nanval;

if ~isfield(cfg,'cntrl_nan'),
    cfg.cntrl_nan = 1;
end

var_names = cfg.vars;

for vIDX=1:size(labels,2),

 		vvals = unique(labels(:,vIDX));
		
		if cfg.cntrl_nan,
			vvals = vvals(vvals~=nanval);
        end
        
	lpre = labels;
	lpre2 = labels;
	
	lpre(:,vIDX) = inf;
	
	ia = {};
	for valIDX=1:length(vvals),
		[ia{valIDX} ib{valIDX} ic{valIDX}] = unique(lpre(labels(:,vIDX)==vvals(valIDX),:),'rows');
	end

	overlap = [];
	for iaIDX=1:length(ia),
		for iaIDX2=(iaIDX+1):length(ia),
			overlap(iaIDX,iaIDX2) = sum(ismember(ia{iaIDX},ia{iaIDX2},'rows'));
		end
    end
    
	if sum(overlap(:))==0,
		excl_var(vIDX) = 1;
		labels(:,vIDX) = inf;
		cs{vIDX} = ia;
	else
		excl_var(vIDX) = 0;
		
		
		ia = {};
		for valIDX=1:length(vvals),
			[ia{valIDX} ib{valIDX} ic{valIDX}] = unique(lpre2(labels(:,vIDX)==vvals(valIDX),:),'rows');
		end
		
		cs{vIDX} = ia;
		
	end
end

for vIDX=1:size(labels,2),
	for vvIDX=1:length(cs{vIDX}),
		cs{vIDX}{vvIDX} = cs{vIDX}{vvIDX}(:,find(sum(isinf(labels))==0));
	end
end
labels = labels(:,find(sum(isinf(labels))==0));
	
stratmat = labels;
[stratcombs, ~,  stratid] = unique(stratmat,'rows');

dm = single([]);
for i=1:size(stratcombs,1),
	dm(:,i) = single(stratid==i);
end

for csIDX=1:length(cs),
	nL = length(cs{csIDX});
	
	cstmp = single([]);
	for lIDX=1:nL,
		cstmp = cat(2,cstmp,ismember(stratcombs,cs{csIDX}{lIDX},'rows'));
	end
	
	cvals = linspace(-1,1,nL);
	cvals = cvals/(sum(abs(cvals))*0.5);
	
	cs_final{csIDX} = cstmp*cvals';
end

cs = cs_final;

end
				

function [params,CCs] = fs_manova_train(data,X,C,varargin),

id_l = find(strcmp(varargin,'lambda'));
if ~isempty(id_l),
	lambda = varargin{id_l+1};
else
	lambda = 0;
end

id_m = find(strcmp(varargin,'m'));
if ~isempty(id_m),
	m = varargin{id_m+1};
else
	m = 1;
end

id_CCs = find(strcmp(varargin,'CCs'));
if ~isempty(id_CCs),
	CCs = varargin{id_CCs+1};
else
	CCs = {};
end

betaDelta = {};

id_B = find(strcmp(varargin,'beta'));
if ~isempty(id_B),
	beta = varargin{id_B+1};
else
	beta = pinv(X)*data;
end
xis = data - X*beta;

if ~iscell(C),
    C = {C};
end

if isempty(CCs),
    for cIDX=1:length(C),
        
        
        CCs{cIDX} = pinv(C{cIDX})*C{cIDX};
        betaDelta{cIDX}      = CCs{cIDX} * beta;
        
    end
else
    for cIDX=1:length(C),
        
        betaDelta{cIDX}      = CCs{cIDX} * beta;
        
    end
end
E = xis'*xis;


E = (1-lambda)*E+lambda*(diag(diag(E)));
iE = pinv(E);

params.betaDelta = betaDelta;
params.iE        = iE;
params.n         = size(data,1); 
params.p         = size(data,2);
params.lambda    = lambda;
params.fE        = params.n-rank(X);
if nargout>1,
    params.CCs = [];
else
    params.CCs = CCs;
end
params.m = m;
params.E = E;
end


function [D,params,params_test,D_st] = fs_manova_test(params,data_test,X_test,varargin),

betaDelta_test = []; cfact = [];


if ~isfield(params,'singleTrial'),
    params.singleTrial = 0;
end

if isfield(params,'beta_test') && ~params.singleTrial,
    beta_test = params.beta_test;
    params_test.beta_test = beta_test;
elseif (~isfield(params,'betaDelta_test') && ~isfield(params,'XX_betaDelta_test')) || params.singleTrial,
    XI = pinv(X_test);
    beta_test = XI*data_test;
    params_test.beta_test = beta_test;
end

if isfield(params,'XX'),
    XX = params.XX;
elseif ~isfield(params,'XX_betaDelta_test'),
    XX = X_test'*X_test;
else
    XX = [];
end


if ~isfield(params,'cfact'),
    params.cfact = [];
end

if isfield(params,'n_used'),
    params_test.n_used = params.n_used;
end

for cIDX=1:length(params.CCs),
    params_test.n_used(cIDX,1) = sum(sum(X_test*params.CCs{cIDX}~=0,2)>0);
    
end

for cIDX=1:length(params.CCs),
    
    if isfield(params,'betaDelta_test') &&  ~isfield(params,'XX_betaDelta_test'),
        betaDelta_test{cIDX} = params.betaDelta_test{cIDX};
        XX_betaDelta_test{cIDX} = XX*betaDelta_test{cIDX};
        
        D(cIDX,1) = trace((params.betaDelta{cIDX}'*((XX*betaDelta_test{cIDX})))*(params.iE));
    elseif ~isfield(params,'XX_betaDelta_test'),
        
        betaDelta_test{cIDX} = params.CCs{cIDX} * beta_test;
        XX_betaDelta_test{cIDX} = XX*betaDelta_test{cIDX};
        
        D(cIDX,1) = trace((params.betaDelta{cIDX}'*((XX_betaDelta_test{cIDX})))*(params.iE));
    else
        D(cIDX,1) = trace((params.betaDelta{cIDX}'*((params.XX_betaDelta_test{cIDX})))*(params.iE));
    end
    
    if params.singleTrial,
        
        bt_all = beta_test;
        bd_all = betaDelta_test{cIDX};
        normbd = [];
        for ccIDX=1:size(bd_all,1),
            normbd(ccIDX) = norm(bd_all(ccIDX,:)).^2;
        end
        for trIDX=1:size(X_test,1),
           
           bt = (XI(:,trIDX))*data_test(trIDX,:);
           XX_t = X_test(trIDX,:)'*X_test(trIDX,:);
           
           idxCnd = find(X_test(trIDX,:));
           bt_c = bt;
           bt_c(:,:) = 0;
           bt_c(idxCnd,:) = bt(idxCnd,:)*params_test.n_used(cIDX);
           
           bd_c   = params.CCs{cIDX}*bt_c;
           
           bd_stp = zeros(size(bd_all));
           for ccIDX=1:size(bd_all,1),
               bd_stp(ccIDX,:) = (dot(bd_all(ccIDX,:),bd_c(ccIDX,:))/normbd(ccIDX))*bd_all(ccIDX,:);
           end
           
           if sum(isnan(bd_stp(:,1)))~=size(bd_stp,1),
               bd_stp(isnan(bd_stp))= 0;
           end
           
           xbd_st = XX*bd_stp;
           
           D_st(cIDX,trIDX) = trace((params.betaDelta{cIDX}'*(squeeze(xbd_st)))*(params.iE));
           D_st(cIDX,trIDX) =  D_st(cIDX,trIDX)*(trace(XX.*XX_t)./trace(XX))*max(sum(params.CCs{cIDX}~=0,2));
          
        end
        
    end
    
end

if isempty(params.cfact),
    cfact =  ((params.fE - params.p - 1)./(params_test.n_used));
else cfact = params.cfact;
end
D = D.*cfact;

idxbad = find((imag(cfact)>0)+(cfact<0));
D(idxbad) = NaN;

if params.singleTrial,
    D_st = D_st.*repmat(cfact,1,size(D_st,2));
    D_st(:,idxbad) = NaN;
    
    for cIDX=1:size(D_st,1),
        D_st(cIDX,:) =  D_st(cIDX,:).*(D(cIDX,:)/(nanmean(D_st(cIDX,:))));
    end
end

params_test.betaDelta_test = betaDelta_test;
params_test.XX_betaDelta_test = XX_betaDelta_test;
params_test.XX = XX;
params_test.cfact = cfact;

end


function [data coeff latent prc_var] = do_pca(condav,data,sid,nComps,coeff,latent)

if nargin<=4,
    if condav,
        dtm = [];
        for cndIDX=1:max(sid),
            dtm(1:size(data,1),:,cndIDX) = mean(data(:,:,sid==cndIDX),3);
        end
        dtm = reshape(dtm,size(dtm,1),[]);
    else
        dtm = reshape(data,size(data,1),[]);
        dtm = dtm(:,randperm(size(dtm,2),min(size(dtm,2),1000000)));
    end
    
    [coeff score latent] = pca(dtm');
end

dt = reshape(data,size(data,1),[]);
ddim = size(data);
if ddim(1)>size(coeff,2),
    ddim(1) = size(coeff,2);
end
dt = dt'*coeff;
dt = reshape(dt',ddim);

dt = dt(1:nComps,:,:);

data = dt;

prc_var = sum(latent(1:nComps))./sum(latent);

end