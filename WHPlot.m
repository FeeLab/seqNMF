function WHPlot(W,H,Data, plotAll, ExtraMatrixToPlot, plotWflat); 
% plots seqNMF factors W and H
% Also plots Data if provided, and reconstruction if data is not provided
% plotAll=1 means plot all data
% plotAll=0 means crop data so it's on the same timescale as W's 
% Emily Mackevicius 2/3/2018
clf
set(gcf, 'color', 'w');
if nargin < 3
    Data = [];
end
if nargin < 4
    plotAll = 1;
end
if nargin<5
    ExtraMatrixToPlot = [];
%     plotAll = 1; % Plotting all the data, and assuming ExtraMatrixToPlot lasts the same duration as Data
end
if nargin<6
    plotWflat = 1;
end

[N,K,L] = size(W); 
[~,T] = size(H);
color_palet = [[0 .6 .3]; [.7 0 .7]; [1 .6 0];  [.1 .3 .9];  [1 .1 .1];  [0 .9 .3]; [.4 .2 .7]; [.7 .2 .1]; [.1 .8 1 ]; [1 .3 .7]; [.2 .8 .2]; [.7 .4 1]; [.9 .6 .4]; [0 .6 1]; [1 .1 .3]]; 
color_palet = repmat(color_palet, ceil(K/size(color_palet,1)),1); 
kColors = color_palet(1:K,:); 
% cleanerW = W.*(W>max(W(:))*.1);
% cleanerH = H.*(H>max(H(:))*.1); 
%% set widths of subplots
m = .1; % margin
ww = .2; % width of W plot
if plotWflat
    wwflat = .05; % width of Wflat plot
else 
    wwflat = 0;
end
hh = .25; % height of H plot
hdata = 1-hh-2*m; 
wdata = 1-ww-wwflat-2*m; 
sep = ceil(L*.5); 

%% crop data, unless plotAll
if plotAll
    indplot = 1:T;
else
    indplot = 2*L+(1:ceil((K*(L+sep))/ww*wdata)); % so that data and W's are on same scale
    indplot(indplot>T) = [];
    TE = size(ExtraMatrixToPlot,2); 
    if TE>0
        ExtraMatrixToPlot(:,[1:floor(indplot(1)/T*TE) ceil(indplot(end)/T*TE):end])=[];
    end
end


%% plot W's
axW = subplot('Position', [m m ww hdata]);cla
hold on
set(gca, 'ColorOrder', kColors); 

WsToPlot = zeros(N,K*(L+sep)); 
XsToPlot = zeros(3,K); 
YsToPlot = [zeros(1,K); N*ones(2,K)]+.5;
for ki = 1:K
    XsToPlot(:,ki) = [(L+sep)*(ki-1) (L+sep)*(ki-1) (L+sep)*(ki-1)+L+1];
    WsToPlot(:,((L+sep)*(ki-1)+1):((L+sep)*(ki-1)+L)) = squeeze(W(:,ki,:));    
end
plot(XsToPlot,YsToPlot, 'linewidth', 2);
clims = [0 prctile(WsToPlot(WsToPlot>0),99)]; % if all W's are empty this line will bug
helper.grayscalePatchPlot(flipud(WsToPlot));%, clims(2)); 



xlim([0 K*(L+sep)]);

if length(ExtraMatrixToPlot)>0
    ylim([0 1.12*N+4])
else
    ylim([0 N+2.5])
end

set(gca, 'ydir', 'normal')
axis off

%% plot data
axIm = subplot('Position', [m+ww m wdata hdata]);cla
if length(Data)~=0
%     clims = [0 prctile(Data(Data>0),99)]; 
    helper.grayscalePatchPlot(flipud(Data(:,indplot)));
%     Im = imagesc(Data(:,indplot), clims);
else
    toplot = helper.reconstruct(W,H(:,indplot));
%     clims = [0 prctile(toplot(:),99.9)]; 
cla
%     for ki = 1:K
%         toplot = helper.reconstruct(W(:,ki,:),H(ki,indplot));
%         helper.grayscalePatchPlot(flipud(toplot), [], kColors(ki,:));
%     end
    helper.grayscalePatchPlot(flipud(toplot));
%     Im = imagesc(toplot,clims);
end
set(gca,'ydir','normal')
hold on; plot([1 1 length(indplot) length(indplot) 1], [0 N N 0 0]+.5, 'k')
xlim([0 length(indplot)+1]); ylim([0 N+2.5])
axis off

if length(ExtraMatrixToPlot)>0
    ExtraMatrixToPlot = ExtraMatrixToPlot - prctile(ExtraMatrixToPlot(:),50); 
    ExtraMatrixToPlot(ExtraMatrixToPlot<=0)=0;
    imagesc(ExtraMatrixToPlot, 'ydata', [N+2.5 1.12*N+2.5], 'xdata', [1 length(indplot)]); 
    cmap  = 1/256*flipud([158,1,66;213,62,79;244,109,67;253,174,97;254,224,139;255,255,191;230,245,152;171,221,164;102,194,165;50,136,189;94,79,162;1 1 1]); % from colorbrewer spectral
    cmap(1,:) = 0; 
    colormap(cmap)
    ylim([0 1.12*N+4])
end

%% plot Wflat (collapse out L dimension of W)
if plotWflat
    axWflat = subplot('Position', [m+ww+wdata m wwflat hdata]);cla
    hold on
    set(gca, 'ColorOrder', kColors); 
    plot(squeeze(sum(W,3)), N:-1:1,'>', 'markersize', 2.5);
    axis tight
    if length(ExtraMatrixToPlot)>0
        ylim([0 1.12*N+4])
    else
        ylim([0 N+2.5])
    end

    xlims = xlim; 
    xlim([xlims(2)*.1 xlims(2)])
    set(gca, 'ydir', 'normal')
    axis off
end
%% plot H's
axH = subplot('Position', [m+ww m+hdata wdata hh]);cla
Hrescaled = repmat(squeeze(sum(sum(W,1),3))',1,T).*H; % rescale by approximate loading
dn = prctile(Hrescaled(:)+eps,100)/2;
for ki = K:-1:1
    Xs = [1 1:length(indplot) length(indplot)]; 
    Ys = [dn*ki (dn*ki + Hrescaled(K-ki+1,indplot)) dn*ki]-dn/2;
    patch(Xs,Ys, kColors(K-ki+1,:), 'edgecolor', kColors(K-ki+1,:))
    hold on
end
ylim([dn/2 dn*K+dn*3]);xlim([0 length(indplot)+1])
axis off
%%
if plotAll
%       linkaxes([axIm axW axWflat], 'y'); 
      linkaxes([axIm axH], 'x');
end
drawnow
