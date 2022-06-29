% run_spikes_stats      run spikes_stats and spikes2spikes for a single directory
%
% CALL                  [ s2s sst ] = run_spikes_stats( filebase )
%
% defaults for optional arguments: 
% clustr        'clc'
% Overwrite     1
% periods       -1: ignore periods during *val* files
% shanknums     []
% logFlag       1
% 
% CALLS                 ParseArgPairs, spikes2spikes, spikes_stats

% 08-jun-12 ES

% revisions
% 16-jun-12 shanknums argument added (otherwise taken from xml)
% 20-aug-13 filebase format
% 23-feb-18 cleanup
% 01-aug-19 (1) call to spikes2spikes and spikes_stats to parameter/value pairs
% 17-aug-19 cleaned up
% 08-oct-19 log file moved to dat directory

function [ s2s, sst ] = run_spikes_stats( filebase, varargin )

s2s = [];
sst = [];

[ clustr, Overwrite, periods, shanknums, logFlag ] = ParseArgPairs( ...
    { 'clustr', 'Overwrite', 'periods', 'shanknums', 'logFlag' } ...
    , { 'clu', 1, -1, [], 1 }...
    , varargin{ : } );

[ pathname, filename, extname ]         = fileparts( filebase );
filename                                = [ filename extname ];
if logFlag
    dn                                  = strfind( filebase, 'dat' );
    if ~isempty( dn )
        logfilename                     = sprintf( '%s/spikes_stats_%s_log.txt', filebase( 1 : ( dn + 2 ) ), filename );
        diary( logfilename )
    end
end

s                                       = dir( sprintf( '%s*%s*', filebase, clustr ) );
if isempty( s )
    fprintf( 1, 'No %s files in %s\n', clustr, pathname );
    return
end

fprintf( 1, 'spikes2spikes for %s...\n', filebase )
s2s                                 = spikes2spikes( filebase, 'clustr', clustr, 'Overwrite', Overwrite, 'periods', periods, 'shanknums', shanknums );
fprintf( 1, 'spikes_stats for %s...\n', filebase )
sst                                 = spikes_stats( filebase, 'shanknums', shanknums, 'clustr', clustr, 'Overwrite', Overwrite, 'periods', periods );

if logFlag && ~isempty( dn )
    diary off
end

return

% EOF
