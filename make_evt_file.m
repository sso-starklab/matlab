% make_evt_file             save *evt file (ndm format - see Hazan et al., 2006, JNM)
%
% call                      rc = make_evt_file( evtfname, mat, spkFs )
%
% gets                      evtfname        full path and name: *.evt.###
%                           mat             [spkFs] matrix of ranges
%                           spkFs           [Hz], 20000
%
% evtfname example
%                          	suf             = 's001'
%                          	evtfname        = sprintf( '%s.evt.%s', filebase, suf );
%
% returns                 	rc              return code from fclose 
%
% calls                     sortranges

% 08-oct-20 ES

function rc = make_evt_file( evtfname, mat, spkFs )

%------------------------------------------------------------------------
% initialize, constants, and arguments
rc                              = -1;
mfname                          = mfilename;

nargs                           = nargin;
if nargs < 2 || isempty( evtfname ) || isempty( mat )
    error( 'missing arguments' )
end
if nargs < 3 || isempty( spkFs )
    spkFs                       = 20000;
end
if size( mat, 2 ) ~= 2 || ~isequal( mat, sortranges( mat ) )
    error( 'mat must be a matrix of ranges' )
end


%------------------------------------------------------------------------
% actual work
evtime                      = mat( :, 1 : 2 ) / spkFs * 1000;               % convert to ms
evlabel                     = ones( size( evtime, 1 ), 1 ) * [ 1 2 ];       % assign labels: 1-onset; 2-offset
fid                         = fopen( evtfname, 'wt' );
fprintf( fid, '%10.2f  %2.0f\n', [ evtime( : )'; evlabel( : )' ] );         % full resolution <=100h of recording...
rc                          = fclose( fid );
fprintf( 1, '%s: Written evt file %s (%d events)\n', mfname, evtfname, size( mat, 1 ) )

return

% EOF
