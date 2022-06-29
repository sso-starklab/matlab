function loopdir
datadir         = '/media/shirly/C22865A128659567/mice/EPS/sst/depth/';
dirname         = [ datadir 'sst/' ];
dirs            = dir( [ dirname '*.sst' ] );
 

nfiles          = length( dirs );
for i           = 1 : nfiles
   fname        = dirs( i ).name;
   sstfname            = [ dirname '.sst' ];
   sstfname     = [ dirname fname ];
   load( sstfname, '-mat' )
%   sst                 = spikes_stats_depth(sstfname , 'graphics', 1,  'Overwrite', 0, 'flipLFP', 1 );
    sst          = spikes_stats_supplement( sst );
   save( sstfname, 'sst' )
end
return


datadir         = '/media/shirly/C22865A128659567/mice/EPS/';
dirname         = [ datadir 's2s/' ];
dirs            = dir( [ dirname '*.s2s' ] );
 

nfiles          = length( dirs );
for i           = 1 : nfiles
   fname        = dirs( i ).name;
   %sstfname            = [ filebase '.sst' ];
   s2sfname     = [ dirname fname ];
   load( s2sfname, '-mat' )
   s2s          = spikes2spikes_supplement( s2s );
   save( s2sfname, 's2s' )
end


return