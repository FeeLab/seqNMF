function p = boxyviolin(Data, binedges, xpos,dx, color)
%%
[N,edges] = histcounts(Data, binedges); 
% below, y means x and x means y :/
Ypts = N; 
Ypts = [0; reshape([Ypts(:) Ypts(:)]',2*length(Ypts),1); 0];
Ypts = Ypts/max(Ypts+eps)*dx/2; 
Xpts = edges; 
Xpts = reshape([Xpts(:) Xpts(:)]',2*length(Xpts),1);
p = patch([xpos + Ypts(:); flipud(xpos - Ypts(:))],[Xpts(:); flipud(Xpts(:))],color, ...
    'linestyle', 'none')
shg