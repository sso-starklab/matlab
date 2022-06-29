% makesrslen        MATLAB version of the bash script of the same name
% 
% [ srslen, fnums, fnames ] = makesrslen( filebase, suffix, wflag )
%
% gets              filebase    full path or par structure
%                   suffix      {'eeg'}
%                   wflag       1       compute and write
%                               0       compute
%                               {-2}    load if exists
%
% does: 
% gives the duration of each file in samples
%
% files:
% input:    *srs file (in case of a merged directory)
%           *xml file (in any case)
% output:   *srslen file (in the parent directory)
%
% calls: /usr/local/bin/makesrslen or get_source_filenums.m

% 09-dec-12 ES

% revisions
% 31-mar-13 load only option (wflag<0)
% 13-may-14 support for merged file only

function [ srslen, fnums, fnames ] = makesrslen( filebase, suffix, wflag )

mfname = upper( mfilename );
nargs = nargin;
if nargs < 1 || isempty( filebase )
    fprintf( 1, '%s: filebase must be supplied\n', mfname );
    return;
end
if nargs < 2 || isempty( suffix )
    suffix = 'eeg';
end
if nargs < 3 || isempty( wflag )
    wflag = -2;
end
if isa( filebase, 'struct' ) && isfield( filebase, 'nBits' )
    par = filebase;
    filebase = par.FileName;
end

srslen = [];
fnums = [];
fnames = [];

[ pathname filename ] = fileparts( filebase );
basedir = fileparts( pathname );
srslenfname = sprintf( '%s.srslen', pathname );

if nargout > 1
    [ fnums fnames ] = get_source_filenums( filebase );
end

if wflag < 0 && exist( srslenfname, 'file' )
    srslen = load( srslenfname );
    if isempty( srslen )
        wflag = 1;
    else
        return
    end
end
    
if isunix && exist( '/usr/local/bin/makesrslen', 'file' ) && strcmp( suffix, 'eeg' ) && wflag && exist( srslenfname, 'file' )
    % use the bash routine
    cd0 = pwd;
    cd( basedir );
    cmd = sprintf( '!/usr/local/bin/makesrslen %s', filename );
    eval( cmd );
    cd( cd0 );
    srslen = load( srslenfname );
elseif exist( pathname, 'dir' )
    % do this locally in MATLAB
    % get the xml file parameters
    if exist( 'par', 'var' ) && isa( par, 'struct' )
        nBits = par.nBits;
        nchans = par.nChannels;
    else
        xmlfname = [ filebase '.xml' ];
        rxml = xmltools( xmlfname );
        rxml = rxml.child(2);
        j = 1;
        while ~strcmp( rxml.child( j ).tag, 'acquisitionSystem' ),
            j = j + 1;
        end
        nBits = str2num( rxml.child( j ).child( 1 ).value);
        nchans = str2num( rxml.child( j ).child( 2 ).value);
    end
    % get the file durations
    if isempty( fnums )
        [ fnums fnames ] = get_source_filenums( filebase );
    end
    nfiles = length( fnums );
    srslen = zeros( nfiles, 1 );
    for i = 1 : length( fnames )
        if length( fnames ) == 1 % a merged dir
            eegfname = [ pathname '/' fnames{ 1 } '.eeg' ];
        else
            eegfname = sprintf( '%s/%s/%s.%s', basedir, fnames{ i }, fnames{ i }, suffix );
        end
        info = dir( eegfname );
        if ~isempty( info ) && isfield( info, 'bytes' )
            srslen( i ) = info.bytes / ( nBits / 8 ) / nchans;
        end
    end
    % save if requested
    if wflag && sum( srslen ) > 0
        fid = fopen( srslenfname, 'w' );
        fprintf( fid, '%d\n', srslen' );
        fclose( fid );
    end
end

return

% EOF

% the bash script is:

% msb138-dl1:dat eranstark$ cat /usr/local/bin/makesrslen 
% #!/bin/bash
% 
% # version 06-dec-12
% 
% echo "making srslen file" 
% 
% size=`ls -l $1.srslen|awk '{print $5}';`; 
% if [ $size -gt 0 ]; then rm $1.srslen; else echo "new file"; fi; 
% #rm $1.srslen
% 
% for i in `cat $1/$1.srs`
% do
% nbits=$(/usr/local/bin/xpathReader $i/$i.xml "//acquisitionSystem/nBits")
% wordlen=$(($nbits/8));
% nchans=$(/usr/local/bin/xpathReader $i/$i.xml "//acquisitionSystem/nChannels")
% nbytes=`ls -l $i/$i.eeg|awk '{print $5}';`
% tmp=$(($nbytes/$nchans));
% nsamples=$(($tmp/$wordlen));
% echo $nsamples>> $1.srslen
% done


% size=`ls -l $1.srslen|awk '{print $5}';`; if [ $size -gt 0 ]; then rm $1.srslen; else echo "new file"; fi; 