function fs_plot_channels(tm,dat),

dat = dat-repmat(mean(dat,2),1,size(dat,2));

dtbins = round(linspace(1,size(dat,2),100));
for binIDX=1:length(dtbins)-1,
    chanrange_act(binIDX) = mean(diff([prctile(dat(:,dtbins(binIDX):dtbins(binIDX+1)),10,2) prctile(dat(:,dtbins(binIDX):dtbins(binIDX+1)),90,2)]'));
end


offset = mean(chanrange_act)*4;

for i=1:size(dat,1),
    hold on;
    plot(tm,dat(i,:)+offset*(i-1));
end

xlim([tm(1) tm(end)]);
ylim([-offset offset*size(dat,1)]);
end