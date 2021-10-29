function vdat = fmcw_burst_snip_signal(vdat, istart, iend)

% creates a function that chops data by istart, iend to prevent false
% signals, a result of voltage hitting voltage cut offs (need to add
% automatic way of doing this
% inspired by Vankova 2018
% EC March 2019

% istart = starting index
% iend = ending index
% vdat = structure from fmcw_load

    vdat.vif = vdat.vif(:,istart:iend); % remove last sample
    disp(['cropped chirps to ' int2str(istart) ':' int2str(iend)])
    vdat.processing = {[mfilename ': cropped chirps to ' int2str(istart) ':' int2str(iend)]};
    vdat = fmcw_derive_parameters(vdat);
    vdat.NSamples = iend-istart;
    vdat.t = vdat.t(istart:iend);
    vdat.f = vdat.f(istart:iend);
    
    vdat = fmcw_burst_make_length_odd(vdat); 
    
end