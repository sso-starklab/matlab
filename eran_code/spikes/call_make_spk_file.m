% call_make_spk_file                generate detrended *spk* files for a given session
%
% call                              rc = call_make_spk_file( filebase )
%
% gets                              filebase
%
% optional arguments (given as name/value pairs):
%
%                                   verbose, detrendflag        see make_spk_file
%                                   Overwrite                   {0} 0: if spk file exists, does not overwrite
%                                                                   1: overwrites in any case
% returns                           rc                          logical (see make_spk_file)
%
% calls                             make_spk_file, LoadXml, LoadRes

% 14-oct-19 ES

function rc = call_make_spk_file( filebase, varargin )

% output initialization
rc                      = [];

% argument handling
nargs = nargin;
if nargs < 1 || isempty( filebase )
    return
end
[ detrendflag, verbose, Overwrite ] = ParseArgPairs(...
    { 'detrendflag', 'verbose', 'Overwrite' } ...
    , { 1, 0, 0 } ...
    , varargin{ : } );

% xml and dat file parameters
par                     = LoadXml( filebase );
nChannels               = par.nChannels;
ng                      = length( par.SpkGrps );
if ng == 0
    fprintf( 1, '%s: no spikes groups in %s.xml\n', upper( mfilename ), filebase )
    return
end
[ pathname, filename ]  = fileparts( filebase );
datfname                = sprintf( '%s/%s.dat', pathname, filename );
if ~exist( datfname, 'file' )
    fprintf( 1, '%s: missing file %s\n', upper( mfilename ), datfname )
    return
end

% go over spike groups and extract spikes
rc                      = zeros( ng, 1 );
for i                   = 1 : ng
    % construct file names
    resfname            = sprintf( '%s/%s.res.%d', pathname, filename, i );
    spkfname            = sprintf( '%s/%s.spk.%d', pathname, filename, i );
    if exist( spkfname, 'file' ) && Overwrite ~= 1
        fprintf( 1, '%s: existing file %s - skipped!!\n', upper( mfilename ), spkfname )
        continue
    end
    % load the res file
    res                 = LoadRes( resfname );
    % prepare arguments for make_spk_file call
    nsamples            = par.SpkGrps( i ).nSamples;
    peaksample          = par.SpkGrps( i ).PeakSample;
    channels            = par.SpkGrps( i ).Channels + 1;
    % call make_spk_file
    rc( i )             = make_spk_file( datfname, spkfname, res, nChannels, channels ...
        , 'nsamples', nsamples, 'peaksample', peaksample ...
        , 'detrendflag', detrendflag, 'verbose', verbose );
end

return

% EOF
