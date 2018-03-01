clf; hold on
K = 10; 
color = .7*[1 1 1]; 
dx = .9; 
binedges = (0:5*K)+.5;
for ki = 1:K
    Data = poissrnd(2*ki,1,200); 
    xpos = ki; 
    p = boxyviolin(Data, binedges, xpos,dx, color)
    Data = randn(1,200)+3*K; 
    xpos = ki; 
    p = boxyviolin(Data, binedges, xpos,dx, [1 .8 .8])
end
shg