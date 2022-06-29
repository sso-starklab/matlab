r_part = [];
kidx_cell = [1;1;0;0;0;0;0;1;1;1;0;1;0;0;1;0;0];
kidx_cell = logical (kidx_cell);
kidx_LED_2 = [0;0;1;1;1;0;0;1;1];
kidx_LED_3 = [0;0;1;0;1;1;0;1;1];
kidx_LED_2 = logical (kidx_LED_2);
kidx_LED_3 = logical (kidx_LED_3);

  [r1_part rc ] = struct_select( r1, kidx_cell,1);

 [r5_part2 rc ] = struct_select( r5, kidx_cell,1);
  [r4_part2 rc ] = struct_select( r4, kidx_cell,1);
   [r3_part2 rc ] = struct_select( r3_part, kidx_LED_3,2);
  [r2_part2 rc ] = struct_select( r2_part, kidx_LED_2,2);
r2_part2.nchans = repmat(length(r2_part2.pAct),1,length(r2_part2.pAct));
r3_part2.nchans = repmat(length(r3_part2.pAct),1,length(r3_part2.pAct));


[ shankclu, fpeth, sgain, pAct ]= seqsort( r4_part2, r5_part2, 'cmp', 'intersect', 'graphics', 1 )

filebase                      = '/media/shirly/C22865A128659567/mice/mC400/dat/mC400_21/mC400_21';
    par                     = LoadXml( filebase );
        amchans                 = get_stimchans( par, [], 'am' );
 [ mmat0, mmat1, am, amTH ] = am2states( filebase, 'channels', amchans, 'nchans', nchans ...
        , 'Fs', Fs, 'TH', amTH, 'graphics', amGraphics, 'Overwrite', amOverwrite ...
        , 'windur', amWindow, 'mindurSEC', mindurSECmov, 'minisiSEC', minisiSECmov, 'amT', amRMS );
    
    
tim_start   = 44*60+22+0.760;
start_idx = tim>tim_start;
end_idx = tim<(tim_start+0.1);
trail_idx = start_idx&end_idx;
sum(trail_idx)

% can only be used after debugstop in segmentBehavior ruf ruf
am_place = '/media/shirly/C22865A128659567/mice/mDS2/dat/mDS2_07/mDS2_07.am';
am0 = fopen(am_place ,'r');
fread(am0,inf,'float32' );

figure;
timss = tim( trail_idx );
timss = timss-min(timss);
mean_am = mean(am);
        plot( timss, am( trail_idx ) );
            set( gca, 'box', 'off', 'tickdir', 'out' )
        xlabel('Time[s]')    
                ylabel( 'Accel [m/s^2]' )
                text(0.1,3.5,sprintf('mean session accel = %0.2g',mean_am))

resize = '-bestfit'; pstr = '-dpdf';
figname = '/media/shirly/Data/_Shirly/Lab/eps/new/acc_fig1';
fig = gcf;

figure( fig );
figi = gcf;
figi.Renderer = 'painters';
pause( 0.2 )
print( fig, pstr, figname, resize )

