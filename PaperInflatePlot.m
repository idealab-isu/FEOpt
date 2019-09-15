%% Initialize
%set(0,'defaulttextinterpreter','latex')
set(0,'DefaultTextFontName','Arial')

homeDir = 'Set this to your working directory';

Vref = [29.95, 20.38, 23.61, 18.73]; %mL
Pref = [910  , 430  , 345  , 675 ]*1e-3; %kPa
%670
Subject = {'D0912','D0917','D1017','D1024'};
SubjectN = {'C1','C2','C3','C4'};

%% Read data

C = cell(4,4);

for i = 1:4
    temp = [homeDir,Subject{i},'/Condo/output/Simulation',Subject{i},'_1_HemoDataInflation.xls'];
    
    InflateData = dlmread(temp,'\t');
    C{i,1} = InflateData(:,6); %V
    C{i,2} = InflateData(:,4); %P
end

%% Generate Klotz
kPa2mmHg = 7.50061561303;

for i = 1:4
    if 0
        Pm = Pref(i); % Moving Klotz Curve
        Vm = max(C{i,1});
        
        C{i,3} = C{i,1};
        
        An = 27.8;
        Bn = 2.76;
        
        V0 = Vm * (0.6 - 0.006*Pm);
        
        V30 = V0 + (Vm-V0) / ((Pm/An) ^ (1/Bn));
        
        beta = log10(Pm/30) / log10(Vm/V30);
        
        alpha = 30 / (V30^beta);
        
        C{i,4} = alpha*C{i,3}.^beta;
    end
    
    if 1
        Pm = Pref(i)*kPa2mmHg; % Fixed Klotz curve endpoint
        Vm = Vref(i);
        C{i,3} = linspace(min(C{i,1}),Vm,length(C{i,1}));
        
        An = 27.8;
        Bn = 2.76;
        
        V0 = Vm * (0.6 - 0.006*Pm);
        
        V30 = V0 + (Vm-V0) / ((Pm/An) ^ (1/Bn));
        
        beta = log10(Pm/30) / log10(Vm/V30);
        
        alpha = 30 / (V30^beta);
        
        C{i,4} = (alpha*C{i,3}.^beta)/kPa2mmHg;
    end
end

%% Plot
close all

for i = 1:4
    plot(C{i,1},C{i,2},'k','LineWidth',2)
    hold on
    plot(C{i,3},C{i,4},'Color',[0.5,0.5,0.5],'LineWidth',0.75)
    %text((max(C{i,1})+max(C{i,3}))/2-1.25,max(C{i,2})+0.05,Subject{i},'FontSize',14)
    text((max(C{i,1})+max(C{i,3}))/2-0.4,max(C{i,2})+0.05,SubjectN{i},'FontSize',14)
end

axis([10,32,0,1.2])

legend({'Simulation','Klotz'},'FontSize',12,'Location','northwest')
legend boxoff 

set(gca,'Visible','off') 
axes('Position',get(gca,'Position'),... 
'XAxisLocation','bottom',... 
'YAxisLocation','left',... 
'Color','none',... 
'XTickLabel',get(gca,'XTickLabel'),... 
'YTickLabel',get(gca,'YTickLabel'),...
'TickDir', 'out',...
'XColor','k','YColor','k');

xlabel('Volume (mL)','FontSize',12)
ylabel('Pressure (kPa)','FontSize',12)

axis([10,32,0,1.2])

%% Save Plot

set(gcf, 'PaperUnits', 'inches');
set(gcf,'PaperSize', [6.25 4.25]);
set(gcf, 'PaperPosition', [0.125 0.125 6 4]);
print(gcf,'-dpdf','-r300','DogInflate.pdf');
