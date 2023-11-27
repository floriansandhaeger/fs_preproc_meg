function [data jumps] = fs_correct_jumps(data),

jumps = [];
idx_meg = find(cellfun(@(x) strcmp(x(1),'M'),data.label));

for chanIDX=1:length(idx_meg),
    disp(['looking for jumps in channel ' num2str(chanIDX)]);
    dtact = data.trial{1}(idx_meg(chanIDX),:);
    dtact_corrected = dtact;
    df = [0 diff(dtact)];
    zd = abs(zscore(df));
    thresh = prctile(zd,99)*20;

    idx_jump = find(zd>thresh);


    for jumpIDX=1:length(idx_jump),
        % only fix jump if this is actually the peak in a local
        % neighbourhood
        interval_act = round([idx_jump(jumpIDX)-data.fsample*0.1 idx_jump(jumpIDX)+data.fsample*0.1]);
        interval_act = setdiff(interval_act(1):interval_act(2),idx_jump(jumpIDX));
        if sum(zd(interval_act)>zd(idx_jump(jumpIDX)))==0,

            % data during jump is removed (+/-10ms)
            % data after jump is shifted
            ns_prepost = round(data.fsample*0.01);
            interval_jump = idx_jump(jumpIDX)+[-ns_prepost:ns_prepost];
            interval_jump = setdiff(interval_jump,[1 length(dtact)]);
            interval_post = (interval_jump(end)+1):length(dtact);
            interval_pre = 1:(interval_jump(1)-1);


            dtact_corrected(interval_jump) = dtact_corrected(interval_pre(end));
            dtact_corrected(interval_post) = dtact_corrected(interval_post) - (dtact_corrected(interval_post(1))-dtact_corrected(interval_pre(end)));

            jumps(end+1,:) = [chanIDX idx_jump(jumpIDX)];
            disp(['corrected jump at sample ' num2str(idx_jump(jumpIDX))]);
        end
    end

    data.trial{1}(idx_meg(chanIDX),:) = dtact_corrected;

end