
function save_print (figname)

resize = '-bestfit'; 
pstr = '-dpdf';
fig = gcf;
resize = '-fillpage';
figure( fig );
figi = gcf;
figi.Renderer = 'painters';
pause( 0.2 )
print( fig, pstr, figname, resize )

return