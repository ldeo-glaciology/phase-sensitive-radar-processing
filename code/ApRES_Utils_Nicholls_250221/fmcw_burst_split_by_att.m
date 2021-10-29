function vdats = fmcw_burst_split_by_att(vdat)

% Splits a fmcw radar burst into (up to) 4 groups of chirps where each
% group has common attenuator settings
%
% Craig Stewart
% 2014/5/2

attSetList = unique(vdat.chirpAtt,'stable');
for ii = 1:length(attSetList)
    vdats(ii) = vdat;
    chirplist = find(vdat.chirpAtt == attSetList(ii)); % we need the actual row numbers here, not logical indexing.
    vdats(ii).processing = [vdat.processing {[mfilename ': keeping chirps with att: ' int2str(real(attSetList(ii))) '+' int2str(imag(attSetList(ii))) ]}]; % ', chirps: ' mat2str(vdat.chirpNum(chirplist))
    vdats(ii) = fmcw_burst_subset(vdats(ii),chirplist);
end