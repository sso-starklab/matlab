% for eran:

%line 233 in postprocess_spikes: rewrite as:
cmd                 = sprintf( '!move %s %s', lfpfname, eegfname );

%line 199 in stim_plot: rewrite as:
figfname = [ figname '.' savetype ];
if ispc
    figfname    = replacetok( figfname, '\', '/' );
end
print( fig, [ figfname '.png' ], [ '-d' savetype ] );


% spikes_stats_plot
% in line 93:
delim                   = strfind( filebase, [ filesep 'dat' filesep ] );
% in line 101:
figname             = [ figdir filesep filename '.sst' ];

% eeg2whl, after line 370, add:
if ispc
    figfname    = replacetok( figfname, '\', '/' );
end

% not clear why plotProbeHFOs runs again
% need to fix error in spikePhase