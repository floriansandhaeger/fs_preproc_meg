function out = ms_convn(a,b);
out = convn(a,b,'same')./convn(ones(size(a)),b,'same');
