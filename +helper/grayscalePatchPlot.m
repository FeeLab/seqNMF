function grayscalePatchPlot(Matrix, scaleby,Color,bounds, alpha)
% bounds where file boundaries are. Default width is 10 time points to the
%        right of boundary in blue. One can modify 10 to smaller or larger numbers
%        based on how long your dataset is.
if nargin<3
    Color = [0 0 0];
end
if nargin<2 || length(scaleby)==0
    scaleby = max(Matrix(:))/3; 
end
if nargin<4
    bounds=[];
end
if nargin<5
    alpha = 1;
end
if isempty(Color);Color = [0 0 0];end
if isempty(scaleby);scaleby = max(Matrix(:))/3;end

Matrix = Matrix./(scaleby);
%Matrix(Matrix<.1) = 0; 
[N,T] = size(Matrix);

for ni = N:-1:1
    Xs = [1 1:T T]; 
    dn = 1;
    Ys = [dn*ni (dn*ni + Matrix(ni,:)) dn*ni]-dn/2;
    patch(Xs,Ys, Color, 'edgecolor', 'none', 'linewidth', .1, 'facealpha', alpha)
    hold on
%     plot([Xs(1) Xs(end)], [Ys(1) Ys(end)],'color', (Color+8*[1 1 1])/9)
end
if ~isempty(bounds)
    for d=1:length(bounds)
        patch([bounds(d),bounds(d)+10,bounds(d)+10,bounds(d)],[0.5,0.5,1+N,1+N], [0 0 1], 'edgecolor', 'none', 'linewidth', .1)
    end
end