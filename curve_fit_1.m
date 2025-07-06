clear,  close all, clc
[Rung, Newt, ddA]=Bfunc; warning('off');
Av0=[.3,1.e-10,1.e-10,.6,.1]'; Ar=1000; a=3; t0=0; ti=0.01; tn=500;
mkr={'o','+','*','.','x','_','|','s','d','^','v','>','<','p','h'};
% -------------------------------------------------------------------------
% id2=[1.0,3.1:.1:3.3, 3.5,3.6, 4.3] ; nid2=length(id2); 
% lgn={'AT','PT','iARD','pARD','WPT','Dz','MRD'}; v=(Ar^2-1)/(Ar^2+1);
% clr2={[0.6,0.3,0],'r','c','m','k','b','g'}; cls=-3.1; dU=@(t) DU(t,1,1);
% for j=1:nid2
%     val=mdl(id2(j));j
%     var={Ar, a, dU, @(t) [],id2(j),cls,val,{2,4}};
%     tic();
%     [Nerr{j},Avn2(j,:)]=Newt([],Av0,var{:});
%     tNR(j)=toc();
%     tic();
%     % [tn(j),Avn1(j,:),t(:,j),Av(:,:,j)]=Rung(t0,ti,tn,Av0R,var{:});
%     [tn1(j),Avn1(j,:),t{j},Av{j}]=Rung(t0,ti,tn,Av0,var{:});tn1(j)
%     tRK(j)=toc();
%     for k=1:length(t{j})
%         dAv{j}(k,:)=ddA([],Av{j}(k,:),v,var{2:end});
%     end
% end
% %
% for m=1:2
%     f=figure(m);clf,f.Color='w'; hold on
%     f.Position=[450,380,675,500];
%     Axx=['\rm a_{' num2str(11*(3-m)) '}'];
%     % title([Axx '-Component'],'FontSize',12);
%     for j=1:nid2
%         mkrj=mkr(mod(j-1,15)+1);
%         % plot(t(:,j),Av(:,3*m-2,j),'Marker',mkrj,'Color',clr2(j,:),...
%         %     'MarkerSize',.5,'LineStyle','-','LineWidth',1);
%         plot(t{j},Av{j}(:,3*m-2),'Color',clr2{j},...
%             'LineStyle','-','LineWidth',1.5);
%     end
%     xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',18);
%     ylabel("$"+Axx+"$",'Interpreter','latex','FontSize',18);
%     xlim([0 200]); legend(lgn{:},'Box','off');
%     set(gca,'GridLineStyle','--','TickDir','both','Box','on',...
%         'FontName','Palatino Linotype','FontSize',14);
%     grid on
%     if m>1
%         ylim([0.5,0.8]); yticks(0.5:.1:0.8);
%     else
%         ylim([0.1,0.4]); yticks(0.1:.1:0.4);
%     end
% end
% %
% f=figure(3);clf,f.Color='w'; hold on; grid on
% f.Position=[450,380,675,500];
% for j=1:nid2
%     lwt=1.5; if j==1, lwt=0.5; end
%     plot(t{j},max(abs(dAv{j}),[],2),'Color',clr2{j},...
%         'LineStyle','-','LineWidth',lwt);
% end
% xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',18);
% % ylabel('$\rm max|\frac{\partial a_{ij}}{\partial t}|$',...
% %     'Interpreter','latex','FontSize',18);
% ylabel("$log(err_{RK4})$",'Interpreter','latex','FontSize',18);
% yscale('log');
% xlim([0 200]); ylim(10.^[-7 -1]); 
% legend(lgn{:},'Box','off','FontSize',14);
% line([0,500],10.^[-6,-6],'Color','k','LineWidth',1.5,...
%     'LineStyle','--','HandleVisibility','off'); 
% set(gca,'GridLineStyle','--','MinorGridLineStyle','none','TickDir','both',...
%     'Box','on','FontName','Palatino Linotype','FontSize',14); 
% %
% err=round(Avn1-Avn2,6)./round(Avn1,6)*100        ;err(isnan(err))=0;
% T=array2table(abs(err(:,[1 4 5],1)));T=varfun(@(x) num2str(x, '%.4f'),T);
% T.Properties.RowNames=lgn;T.Properties.VariableNames={'A_11','A_22','A_13'};
% %
% f=figure(4); clf,f.Color='w'; hold on
% f.Position=[450,380,675,500];
% for j=1:nid2
%     plot(Nerr{j},'Color',clr2{j},'Marker','x','MarkerSize',5,...
%         'LineStyle','-.','LineWidth',1.5);
% end
% line([1,8],10.^[-6,-6],'Color','k','LineWidth',1.,...
%     'LineStyle','--','HandleVisibility','off'); 
% xlabel("$Iteration$",'Interpreter','latex'); 
% ylabel("$log(err_{NR})$",'Interpreter','latex');
% yscale('log'); legend(lgn{:},'Box','off','Location','southwest'); grid on
% set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
%     'TickDir','both','Box','on'); 
% T2=array2table([tNR(:),tRK(:)]'); T2.Properties.VariableNames=lgn;
% T2.Properties.RowNames=["NR","RK45"];
% -------------------------------------------------------------------------
% cls={};
% cls{1}=[1.1,1.2,2.1,2.2,2.3,2.4,2.5,2.6]; v=(Ar^2-1)/(Ar^2+1);
% dsn{1}=["HYB_1","HYB_2","ISO","LIN","QDR","SF2","HL1","HL2"];
% 
% cls{2}=[3.0,-1.1,-2.1,-2.3,-2.4,-3.1,-3.2,-3.3,-3.4,-4.1,-4.2,-4.3]; 
% dsn{2}=["IBOF","ORS","ORT","NAT_1","ORW","NAT_2",...
%         "WTZ","LAR32","ORW3","VST","FFLAR4","LAR4"];
% G=1; dU=@(t) DU(t,G,1);ncls(1)=length(cls{1});ncls(2)=length(cls{2});
% %
% Nerr=[];Avn2=[];Avn1=[]; Av=[]; t=[]; 
% clr2{1}=[255,0,  0;0  ,255,0;0  ,0,255;0,255,255;
%          255,0,255;255,165,0;128,0,128;0,  0,  0]/255;
% clr2{2}=[
% 255,  0,0;0  ,255,  0;0  ,0  ,255;255,255, 0;0  ,255,255;255,0,255;
% 255,165,0;128,0  ,128;255,192,203;139, 69,19;128,128,128;  0,0,  0]/255;
% %
% for i=1:2
% nclsi=ncls(i); 
% for j=1:nclsi
%     var={Ar, a, dU, @(t) [], 1.,cls{i}(j),{1., .01},{2,4}};
%     tic()
%     [Nerr{j,i}, Avn2(j,:,i)]=Newt([], Av0,var{:});
%     tNR(j,i)=toc();
%     tic()
%     [tn2(j,i),Avn1(j,:,i),t{i,j},Av{i,j}]=Rung(t0,ti,tn,Av0,var{:});tn2(j,i)
%     tRK(j,i)=toc();
%     for k=1:length(t{i,j})
%         dAv{i,j}(k,:)=ddA([],Av{i,j}(k,:),v,var{2:end});
%     end
% end
% end
% %
% for i=1:2
%     for j=1:2
%         f=figure(j+2*(i-1)+4); clf,f.Color='w'; hold on
%         f.Position=[450,360,600,540];
%         Axx=['\rm a_{' num2str(10+j) '}'];
%         jlen(i)=0;
%         for k=1:ncls(i)
%             mkrk=mkr(mod(k-1,15)+1);
%             plot(G*t{i,k},Av{i,k}(:,j+3),'Marker',mkrk,'Color',clr2{i}(k,:),...
%                 'MarkerSize',1,'LineWidth',.5);
%             jlen(i)=max(jlen(i),max(G*t{i,k}));
%         end
%         xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',18);
%         ylabel("$"+Axx+"$",'Interpreter','latex','FontSize',18);
%         xlim([0 1.05*jlen(i)]);  ylim([0,1]-[.2,.7]*(j-1));
%         grid on; set(gca,'Box','on','Tickdir','both',...
%             'FontName','Palatino Linotype','FontSize',14);
%         legend(dsn{i}(:),'Location','southoutside','NumColumns',5,...
%             'Orientation','horizontal','FontSize',14,'Box','off');
%     end
% end
% %
% for i=1:2
%     f=figure(i+8); clf,f.Color='w'; f.Position=[450,360,600,540]; hold on
%     for k=1:ncls(i)
%         mkrk=mkr(mod(k-1,15)+1);
%         plot(G*t{i,k},max(abs(dAv{i,k}),[],2),'Color',clr2{i}(k,:),...
%             'LineWidth',1.);
%     end
%     line([0,1.05*jlen(i)],10.^[-6,-6],'Color','k','LineWidth',1,...
%         'LineStyle','--','HandleVisibility','off'); 
%     xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',18); 
%     % ylabel('$\rm max|\frac{\partial a_{ij}}{\partial t}|$',...
%     %     'Interpreter','latex','FontSize',18); 
%     ylabel("$log(err_{RK45})$",'Interpreter','latex','FontSize',18); 
%     yscale('log'); ylim(10.^[-7 1-i]); yticks(logspace(-7,1-i,9-i));
%     xlim([0, 1.05*jlen(i)]);
%     set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
%         'TickDir','both','Box','on','FontSize',14,'FontName',...
%         'Palatino Linotype'); grid on; 
%     legend(dsn{i}(:),'Box','off','Location','southoutside','Orientation',...
%    'horizontal','NumColumns',4,'FontSize',14);
% end
% %
% err=round(Avn1-Avn2,5)./round(Avn1,5)*100; 
% err(isnan(err))=0;err=err(:,[1 4 5],:);
% for j=1:2
%     Tj=array2table(abs(err(1:ncls(j),:,j)));
%     T{j}=varfun(@(x) num2str(x, '%.4f'),Tj);
%     T{j}.Properties.VariableNames={'A_11','A_22','A_12'};
%     T{j}.Properties.RowNames=dsn{j};
% end
% for i=1:2
%     f=figure(i+10); clf, f.Color='w'; hold on
%     f.Position=[450,380,675,500];
%     nclsi=ncls(i); ilen=0;
%     for j=1:nclsi
%         plot(Nerr{j,i},'Color',clr2{i}(j,:),'Marker','x','MarkerSize',5,...
%             'LineStyle','-.','LineWidth',1.);
%         ilen=max(ilen,length(Nerr{j,i}));
%     end
%     line([1,ilen],10.^[-6,-6],'Color','k','LineWidth',1.,...
%         'LineStyle','--','HandleVisibility','off');
%     xlabel("$Iteration$"    ,'Interpreter','latex','FontSize',18);
%     ylabel("$log(err_{NR})$",'Interpreter','latex','FontSize',18); 
%     yscale('log');  grid on
%     set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
%         'TickDir','both','Box','on','FontName','Palatino Linotype',...
%         'FontSize',14);
%     legend(dsn{i}(:),'Location','southoutside','NumColumns',5,...
%         'Orientation','horizontal','FontSize',14,'Box','off');
% end
%% ------------------------------------------------------------------------
function f=mdl(flg)
flg1=floor(flg); flg2=round(10*(flg-flg1));
switch flg1
    case 1 % FT - SRF
        f={ 1 .0311};
    case 2 % RSC
        f={.2 .0311};
end
switch flg1
    case {3 4 6 7}
        switch flg2
            case 0 % IRD
                f={.0311}; % 'CI'
            case 1 % ARD
                f={[1.924,58.390,400.,.1168,0.]*1e-4}; % 'bta 1-5'
            case {2 4} % iARD
                f={.0562 .9977};     % ['CI' 'CM'] ;
            case 3 % pARD
                switch flg1
                    case 4 % mARD
                        %D=zeros(3); D([1 5 9])=[1.,.8,.15];f={.04796 D}; % ['CI' 'D' ] ;
                        D=zeros(3); D([1 5 9])=[1.,.7946,.012];f={.0198 D}; % ['CI' 'D' ] ;
                    otherwise % pARD
                        c=.9868; D=[1,c,1-c]; D=diag(D,0);
                        f={.0169 D};      % ['CI' 'D' ] ;
                end
            case 5 % WPT
                f={.0504 .9950}; % ['CI', 'w']
            case 6 % Dz
                n=[0,0,1];
                f={.0258, .0051, n}; % ['CI', 'Dz', 'n']
        end
    case 5 % NEM
        f={.0311; 1};  % ['CI, U0']
end
switch flg1
    case 6 % (p)ARD-RSC
        f=[{1/30} f]  ; %['kpa' ARD-vars]
    case 7 % X-RPR
        f=[{1-1/30, .01} f]; %['lfa' 'bta', X-vars]
end
end
%% ------------------------------------------------------------------------
function f=DU(t,flg,varargin)
f=zeros(3);
switch flg
    case 1 % Simple shear
        G=varargin{1};  f(6)=G;
    case 2 % Shearing/stretching
        [E,G]=varargin{1:2}; 
        f([1,5,6,9])=[2*E,-E,G,-E];
    case 3 % Uniaxial
        E=varargin{1}; 
        f([1,5,9])=[-E,-E,2*E];
    case 4 % Biaxial
        E=varargin{1};
        f([1,5,9])=[-2*E, E, E];
    case 5 % Shearing/planar stretching
        [E,G]=varargin{1:2};
        f([3,5,9])=[G,E,-E];

end
end