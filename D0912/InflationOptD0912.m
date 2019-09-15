%% Initialize Results Log

if exist('results.txt','file')
    copyfile('results.txt','results2.txt')
end

dlmwrite('results.txt',zeros(1,7),'delimiter',',','precision',9);

%% Parameters and Options

x0 = [ 8   ,  1.5 , 15   , 15   ]; %[15, 1  , 25, 0.5 ] 
%     b    , a    , bf   , af
lb = [ 0.01,  0.01,  0.01,  0.01]; %[ 5, 0.5, 10, 0.25]
ub = [50   , 50   , 50   , 50   ]; %[25, 1.5, 40, 1   ]

%% Optimization

tic % Pattern Search
options = optimoptions('patternsearch','Display','iter','UseCompletePoll',true,'UseVectorized',true,'FunctionTolerance',1e-2,'StepTolerance',1e-2,'MeshTolerance',1e-3);
[Xmin,Smin] = patternsearch(@(x)InflateD0912vec(x),x0,[],[],[],[],lb,ub,[],options);
totalTime = toc;

%tic % PSO
%options = optimoptions('particleswarm','Display','iter','UseVectorized',true,'MeshTolerance',5e-4);
%[Xmin,Smin] = particleswarm(@(x)InflateD0912vec(x),size(x0,2),lb,ub,options);
%totalTime = toc;

save('finalD0912.mat')

InflateD0912vec(Xmin,1)

csvwrite('stopFlag.txt',1);
