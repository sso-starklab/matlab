

datadir                 = '/media/shirly/C22865A128659567/mice/EPS/';
filebase                = '';
suffix                  = si.filebase{ 1 };
ctwildcard              = sprintf( '%sdc/%s.*_c%s.celltypeClassification', datadir, suffix, ilevel );


%mDL5
    L = load ('/odin/Shirly/dc/fromAmir/mDL5_32.T33_cB.celltypeClassification', '-mat');
    temp {1}= 'pv::chr2'; % T33 T36 BLUE

    L = load ('/odin/Shirly/dc/fromAmir/mDL5_32.T37T38_cB.celltypeClassification', '-mat');
    temp {1}= 'pv::jaws'; % 34 35 T37 T38 RED

    L.s.opsinType=repmat(temp,32,1)
    z=L;


    save ('/odin/Shirly/dc/mDL5_32.T33_cB.celltypeClassification', 'z')

    save ('/odin/Shirly/dc/mDL5_32.T37T38_cB.celltypeClassification', 'z')

    save ('/odin/Shirly/dc/mDL5_27.T34T35_cB.celltypeClassification', 'z')

%mP23

    L = load ('/odin/Shirly/dc/fromAmir/mP23_29.T33T36_cB.celltypeClassification', '-mat');
    temp {1}= 'pv::chr2'; % T33 T36 BLUE

    L = load ('/odin/Shirly/dc/fromAmir/mP23_29.T37T38_cB.celltypeClassification', '-mat');
    L = load ('/odin/Shirly/dc/fromAmir/mP23_29.T34T35_cB.celltypeClassification', '-mat');

    temp {1}= 'pv::jaws'; % 34 35 T37 T38 RED

    L.s.opsinType=repmat(temp,28,1)
    z=L;


    save ('/media/shirly/C22865A128659567/mice/EPS/dc/mDS1_15.T65T66_cB.celltypeClassification', 'z')

    save ('/odin/Shirly/dc/mP23_29.T37T38_cB.celltypeClassification', 'z')

    save ('/odin/Shirly/dc/mP23_29.T34T35_cB.celltypeClassification', 'z')


% for nested structures
    ctfullname = '/media/shirly/C22865A128659567/mice/EPS/dc/mDL5_06.T33_cB.celltypeClassification';
    L = load( ctfullname, '-mat' );
    s = L.z.s; 
    stats = L.z.stats; 
    save( ctfullname, 's', 'stats' )

