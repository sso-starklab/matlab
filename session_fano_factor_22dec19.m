% session_fano_factor       organize session data and compute Fano factors
%
% call                      [ fanomat, s ] = session_fano_factor( filebase, varargin )
%
% gets                      filebase
%
% optional arguments (given as name/value pairs):
% 
%                           brain_state         {'SWS'}
%                           maxwsMS             {1024}
%                           minwsMS             {1}
%
% does                      (1) determine periods during the requested brain state, without stimuli 
%                           (2) loads spikes of all well-isolated units
%                           (3) keep relevant spikes and re-time to integer multiples of the largest requested window size
%                           (4) call a computational routine (calc_fano_factor)
%                           (5) organize in a structure and plot
% 
% returns                   fanomat             matrix of nclu x nwindows (see calc_fano_factor)
%                           s                   structure with additional details 
%
% calls                     LoadXml                                 (blab)
%                           get_states, LoadStims                   (formats)
%                           ParseArgPairs                           (general)
%                           lin2log                                 (graph)
%                           inranges, resampleranges, setdiffranges (sets)
%                           calc_fano_factor, load_spikes           (spikes)

% 15-dec-19 SSo + ES

% revisions
% 22-dec-19 (1) identified (but not solved yet) issues in steps 6 and 7
%           (2) updated call to calc_fano_factor (step 8)
%           (3) added a conceptual step 9

function [ fanomat, s ] = session_fano_factor( filebase, varargin )

% constants
colors              = [ 0 0 0.7; 1 0 0 ];

% initialize output
fanomat             = [];

% handle arguments
nargs               = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ brain_state, minwsMS, maxwsMS ...
    , padBufferSEC, shanknums, ilevel ...
    , graphics ] = ParseArgPairs(...
    { 'brain_state', 'minwsMS', 'maxwsMS' ...
    , 'padBufferSEC', 'shanknums', 'ilevel' ...
    , 'graphics', }...
    , { 'SWS' 1, 1024 ...
    , [ -0.01 0.01 ], [], 'B' ...
    , 1 } ...
    , varargin{ : } );

% get session parameters
par                 = LoadXml( filebase );
spkFs               = par.SampleRate;
eegFs               = par.lfpSampleRate;

% 1. determine vector of integer-sized, log-spaced window sizes 
wsMS                = [ minwsMS maxwsMS ];
ws                  = wsMS * spkFs / 1000;         % convert requested window to samples
ws                  = 2.^floor( log2( ws ) );       % floor maximal window size to an integer power of 2
wsvec               = 2 .^ ( log2( ws( 1 ) ) : 1 : log2( ws( 2 ) ) );
nwindows            = length( wsvec );

% 2. determine periods of the requested brain state
periods             = get_states( filebase, brain_state );              % get the periods for the brain_state
periods             = resampleranges( periods, spkFs, eegFs );          % resample ranges to spkFs

% 3. clean periods 
vals                = LoadStims( filebase );                            % determine stimulus times
stims               = vals( :, 1 : 2 );
pads                = padBufferSEC * spkFs;
blackout            = stims + ones( size( stims, 1 ), 1 ) * pads;       % pad-buffer stimulus times 
periods             = setdiffranges( periods, blackout );               % exclude stimuli from the periods

% 4. clip all periods to integer multiples of the maxws
maxws               = wsvec( nwindows );
durs                = diff( periods, [], 2 ) + 1;
ridx                = durs < maxws;                                     % remove short periods
periods( ridx, : )  = [];
durs( ridx, : )     = [];
periods             = [ periods( :, 1 ) periods( :, 1 ) + floor( durs / maxws ) * maxws - 1 ]; % actually clip
durs                = diff( periods, [], 2 ) + 1;
nperiods            = size( periods, 1 );
offset              = [ 0; cumsum( durs( 1 : nperiods - 1 ) / maxws ) ] * maxws;    % offset in samples for each spike in each "block"

% 5. load spikes
spk                 = load_spikes( filebase, shanknums, ilevel );
uclu                = unique( spk.clu );
nclu                = length( uclu );

% 6. keep only the spikes during the clipped periods
[ idx, midx, res ]  = inranges( spk.res, periods, 1 );                  % remove time of first sample in each period
clu                 = spk.clu( idx );

% 7. re-time the res to be continuous
reshat              = NaN( size( res ) );
for i               = 1 : nperiods
    pidx            = midx == i;
    reshat( pidx )  = res( pidx ) + offset( i );
end

% 8. call a computation routine:
%[ fanomat, myu, sd2, nwindows, nspikes ] = calc_fano_factor( spk.res, spk.clu, 'winsizes', wsvec );
[ mat1, mat2, mat3, nwindows, vec ] = calc_fano_factor( reshat, clu, 'winsizes', wsvec );

% 9. project fanomat from calc_fano_factor to fanomat with the full nclu
midx                = ismember( uclu, unique( clu ) );
fanomat             = NaN( nclu, length( wsvec ) );
myu                 = fanomat;
sd2                 = fanomat;
fanomat( midx, : )  = mat1;
myu( midx, : )      = mat2;
sd2( midx, : )      = mat3;
nspikes             = NaN( nclu, 1 );
nspikes( midx )     = vec;

% 10. organize detailed output
s.shankclu          = spk.shankclu;
s.fanomat           = fanomat;
s.myu               = myu;
s.sd2               = sd2;
s.nspikes           = nspikes;
s.nwindows          = nwindows;
s.winsizes          = wsvec;

% 11. graphical summary
if graphics
    figure
    for cti         = 1 : 2
        ct          = cti - 1;
        idx         = s.shankclu( :, 3 ) == ct;
        subplot( 2, 2, cti )
        fanomatI    = s.fanomat( idx, : );
        ph          = plot( log2( s.winsizes ), fanomatI', '.-b' );
        set( ph, 'color', colors( cti, : ) )
        lin2log( 'x', 2, 1 )
        set( gca, 'tickdir', 'out', 'box', 'off' ),
        xlabel( 'Window size [samples]' )
        ylabel( 'Fano Factor' )
        ylims = [ min( s.fanomat( : ) ) * 0.9 max( s.fanomat( : ) ) * 1.1 ];
        set( gca, 'yscale', 'log', 'ylim', ylims )
    end
end

return

% EOF

filebase = 'G:\mice\mK01\mK01_10';

savepath = 'G:\mice\EPS\fanofactor';
[ ~, filename ] = fileparts( filebase );
[ fanomat, s ] = session_fano_factor( filebase );
savename = [ savepath '\' filename '_fanofactors.mat' ];
save( savename, 's', 'fanomat' );
