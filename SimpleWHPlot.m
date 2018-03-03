function SimpleWHPlot(W,H,Data, plotAll); 
% plots seqNMF factors W and H
% Also plots Data if provided, and reconstruction if data is not provided
% plotAll=1 means plot all data
% plotAll=0 means crop data so it's on the same timescale as W's 
% Emily Mackevicius 2/3/2018

clf

% set(gcf, 'color', 'w');
if nargin < 4
    plotAll = 1;
end
if nargin < 3 || size(Data,2)==0
    plotData = 0;
else
    plotData = 1; 
end

[N,K,L] = size(W); 
[~,T] = size(H);
color_palet = [[0 .6 .3]; [.7 0 .7]; [1 .6 0];  [.1 .3 .9];  [1 .1 .1];  [0 .9 .3]; [.4 .2 .7]; [.7 .2 .1]; [.1 .8 1 ]; [1 .3 .7]; [.2 .8 .2]; [.7 .4 1]; [.9 .6 .4]; [0 .6 1]; [1 .1 .3]]; 
color_palet = repmat(color_palet, ceil(K/size(color_palet,1)),1); 
kColors = color_palet(1:K,:); 
%% set widths of subplots
m = .1; % margin
ww = .2; % width of W plot
wwflat = .05; % width of Wflat plot
hh = .25; % height of H plot
hdata = 1-hh-2*m; 
wdata = 1-ww-wwflat-2*m; 
sep = ceil(L*.1); 

%% crop data, unless plotAll
if plotAll
    indplot = 1:T;
else
    indplot = 2*L+(1:ceil((K*(L+sep))/ww*wdata)); % so that data and W's are on same scale
    indplot(indplot>T) = [];
end


%% plot W's
axW = subplot('Position', [m m ww hdata]);
hold on
set(gca, 'ColorOrder', kColors); 

WsToPlot = zeros(N,K*(L+sep)); 
XsToPlot = zeros(3,K); 
YsToPlot = [N*ones(1,K); zeros(2,K)]+.5;
for ki = 1:K
    WsToPlot(:,((L+sep)*(ki-1)+1):((L+sep)*(ki-1)+L)) = squeeze(W(:,ki,:));
    XsToPlot(:,ki) = [(L+sep)*(ki-1)+1 (L+sep)*(ki-1)+1 (L+sep)*(ki-1)+L];
end
clims = [0 prctile(WsToPlot(WsToPlot>0),99)]; % if all W's are empty this line will bug
imagesc(WsToPlot, clims); 
% cmap =
% 1/256*flipud([158,1,66;213,62,79;244,109,67;253,174,97;254,224,139;255,255,191;230,245,152;171,221,164;102,194,165;50,136,189;94,79,162;1 1 1]); % colorbrewer
cmap = flipud(gray); 
colormap(cmap)
plot(XsToPlot,YsToPlot);
xlim([1 K*(L+sep)]);ylim([0 N+.1]+.5)
set(gca, 'ydir', 'reverse')
axis off

%% plot data
axIm = subplot('Position', [m+ww m wdata hdata]);
if plotData
    clims = [0 prctile(Data(Data>0),99)]; 
    Im = imagesc(Data(:,indplot), clims);
else
    toplot = helper.reconstruct(W,H(:,indplot));
    clims = [0 prctile(toplot(:),99.9)]; 
    Im = imagesc(toplot,clims);
end
set(gca,'ydir','reverse')
hold on; plot([0 0 length(indplot)+1], [N 0 0]+.5, 'k')
xlim([0 length(indplot)+1]);ylim([0 N+.1]+.5)
axis off
%% plot Wflat (collapse out L dimension of W)
axWflat = subplot('Position', [m+ww+wdata m wwflat hdata]);
hold on
set(gca, 'ColorOrder', kColors); 
plot(squeeze(sum(W,3)), 1:N,'>', 'markersize', 2.5);
ylim([0 N+.1]+.5)
axis tight
xlims = xlim; 
xlim([xlims(2)*.1 xlims(2)])
set(gca, 'ydir', 'reverse')
axis off
%% plot H's
axH = subplot('Position', [m+ww m+hdata wdata hh]);
Hrescaled = repmat(squeeze(sum(sum(W,1),3))',1,T).*H; % rescale by approximate loading
dn = prctile(Hrescaled(:),100)/2; 
for ki = K:-1:1
    Xs = [1 1:length(indplot) length(indplot)]; 
    Ys = [dn*ki (dn*ki + Hrescaled(K-ki+1,indplot)) dn*ki]-dn/2;
    patch(Xs,Ys, kColors(K-ki+1,:), 'edgecolor', kColors(K-ki+1,:))
    hold on
end
ylim([0 dn*K+dn*3]);xlim([0 length(indplot)+1])
axis off
%%
if plotAll
      linkaxes([axIm axW axWflat], 'y'); linkaxes([axIm axH], 'x');
end

