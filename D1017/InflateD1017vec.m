function [y] = InflateD1017vec(x,~)

subject = 'D1017';
Pref = 345; % Pressure
Vref = 23.61; % Volume

homeDir = []; % Set to home directory
contDir = []; % Set to coninuity directory

%%%%%

n = size(x,1) % Number of input points
y = Inf*ones(n,1); % Output
R = zeros(n,3); % Full results
status = zeros(n,1); % Used to track each evaluation's iteration progress % -1 = done, 0 = unprocessed, 1 = running

csvwrite('numWorkers.txt',n);

%% Check for Repeats

if nargin == 1 % Don't ignore repeats
    for i = 1:n
        Ri = repeat(x(i,:)); %Check for repeat [runs,Vm,RMS]
        
        if ~isempty(Ri) % Repeat configuration
            disp(['repeat',num2str(i)])
            
            R(i,:) = Ri;
            y(i) = Ri(3); % RMS
            status(i) = -1; % Don't run this evaluation
        end
    end
end

%% Evaluate

firstRun = 1;

while max(status) >= 0 % Still configurations evaluating
    for i = 1:n
        if status(i) >= 0 % Configuration 'i' not complete
            
            if status(i) == 0 % Not started
                copyfile([subject,'_0_ED_DT_Inflation.cont6'],[subject,'_',num2str(i),'_ED_DT_Inflation.cont6']); % Copy inflation model
                copyfile([subject,'_0_ED_DT_Deflation.cont6'],[subject,'_',num2str(i),'_ED_DT_Deflation.cont6']); % Copy deflation model
                copyfile([subject,'_0_Nodes_ED_DT.xls'],[subject,'_',num2str(i),'_Nodes_ED_DT.xls']); % Copy nodes
                
                dlmwrite(['inputs',num2str(i),'.txt'],x(i,:)','delimiter',',','precision',9);
                
                status(i) = 1; % Configuration has been started
                
                disp(['worker',num2str(i),' ready'])
                
            else % Running
                state = 1;
                try
                    state = csvread(['state',num2str(i),'.txt']);
                catch
                end
                
                if state < 0 % Done
                    %disp(['worker',num2str(i),' complete'])
                    
                    delete([subject,'_',num2str(i),'_ED_DT_Inflation.cont6']); % Clear for next run
                    delete([subject,'_',num2str(i),'_ED_DT_Deflation.cont6']);
                    delete([subject,'_',num2str(i),'_Nodes_ED_DT.xls']);
                    delete(['inputs',num2str(i),'.txt']);
		    delete(['state',num2str(i),'.txt']);
                    
                    if exist([homeDir,subject,'/Inflation_',num2str(i),'_',num2str(-state),'/StressStrain',num2str(Pref),'.xls'], 'file') && exist([contDir,'.continuity/working/Simulation',subject,'_',num2str(i),'_HemoDataInflation.xls'], 'file') % Evaluation succesful
                        
                        Ri = process(i,subject,Vref); % Results processing
                        
                        if isempty(Ri) % Not converged
                            R(i,:) = [state,10*Vref,Inf];
                            y(i) = 1e6;
                            
                            status(i) = -1;
                            
                        else % Converged
                            Ri(1) = status(i);
                            R(i,:) = Ri;
                            y(i) = Ri(3);
                            
                            status(i) = -1;
                        end
                    else
                        disp(['failed',num2str(i)])
                        
                        R(i,:) = [state,0,Inf];
                        y(i) = 1e6;
                        
                        status(i) = -1;
                    end
                end
            end
        end
    end
    
    if firstRun
        csvwrite('readyFlag.txt',1); % All files are prepared
        firstRun = 0;
        disp('all workers ready')
    end
end

%delete(['readyFlag.txt']);

logger(x,R);


%% Check for Repeat
function Ri = repeat(x)

Ri = [];

results = [];

while isempty(results) % Try until the file can be read
    try
        results = csvread('results.txt'); % Read file
    catch
        pause(0.01);
    end
end

for i = 1:size(results,1) % Check results for matching parameters
    if sqrt(sum( (results(i,1:4)-x(1:4)).^2 )) < 0.000001
        Ri = results(i,5:7);
        Ri(1) = -Ri(1);
        break
    end
end

%% Process Results
function R = process(i,subject,Vref)

InflateData = dlmread([contDir,'.continuity/working/Simulation',subject,'_',num2str(i),'_HemoDataInflation.xls'],'\t');

InflateP = InflateData(:,4); %kPa
InflateV = InflateData(:,6); %mL

kPa2mmHg = 7.50061561303;

InflateP = InflateP*kPa2mmHg; %mmHg

Pm = InflateP(end); %mmHg
Vm = InflateV(end); %mL

if abs(Vm - Vref)/Vref > 0.03 % Volume not converged
    disp(['Vm',num2str(i),': ',num2str(Vm),', diverged'])
    R = [];
else % Volume converged
    disp(['Vm',num2str(i),': ',num2str(Vm),', converged'])
    
    An = 27.8;
    Bn = 2.76;
    
    V0 = Vm * (0.6 - 0.006*Pm);
    
    V30 = V0 + (Vm-V0) / ((Pm/An) ^ (1/Bn));
    
    beta = log10(Pm/30) / log10(Vm/V30);
    
    alpha = 30 / (V30^beta);
    
    KlotzP = alpha*InflateV.^beta;
    
    RMS = sqrt(sum((KlotzP-InflateP).^2));
    
    R = [0.5,Vm,RMS];
end

%% Log Result
function logger(x,R)

results = [];

while isempty(results) % Try until the file can be read
    try
        results = csvread('results.txt'); % Read file
    catch
        pause(0.01);
    end
end

if results(1,1) == 0
    results = [x,R];
else
    results = [results;x,R];
end

try
    dlmwrite('results.txt',results,'delimiter',',','precision',9);
catch
end
