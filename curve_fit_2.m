clear,  close all, clc
T0=0; h=1.; Tn=6000; [Rung, Newt, ddA]=Bfunc;
Ar=1000; a=3;CI=.01;CLS=-3.1; v=(Ar^2-1)/(Ar^2+1);
% -------------------------------------------------------------------------
% Av0=[.3,1.e-10,1.e-10,.6,.1]';
% E=.1;G=10*E; dU=zeros(3,3,2); dU(6)=G; dU(9+[1 6 9])=[E G -E];
% id=[2 1 1 7.0] ;nid=length(id); kp={{.1}, {1.}, {.1},{.9, .0}};
% for i=1:2
%     for j=1:nid
%         k=j+(i-1)*nid;
%         var={Ar, a,@(t) dU(:,:,i),@(t) [],id(j),CLS,[kp{j}, {CI}],{2,4}};
%         tic();  [Nerr{j,i},Avn2(j,:,i)]=Newt([],Av0,var{:}); tNR(j)=toc();
%         tic();
%         [tn(j,i),Avn1(j,:,i),t{i,j},Av{i,j}]=Rung(T0,h,Tn,Av0,var{:});k
%         tRK(j)=toc();
%         for rk=1:length(t{i,j})
%             dAv{i,j}(rk,:)=ddA([],Av{i,j}(rk,:),v,var{2:end});
%         end
%     end
% end
% %
% clr={'k', 'g','r','c'}; nms={'RSC', "AT ", 'SRF','RPR'};
% lst={'--', '-', '-.'}; pk=[1 4 5];
% for i=1:2
%     f=figure(i); clf; f.Color='w'; grid on; hold on
%     f.Position=[450,380,675,500];
%     for j=1:nid
%         for k=1:3
%             Axx="$\rm "+nms{j}+"-a_{"+string(5*(pk(k)-k)+10+k)+"}$";
%             p=plot(G*t{i,j},Av{i,j}(:,pk(k)),[clr{j},lst{k}],'DisplayName',Axx,...
%                 'LineWidth',1.5);
%         end
%     end
%     xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',18);
%     ylabel("$\rm a_{ij}$",'Interpreter','latex','FontSize',18);
%     xlim([0 750]);ylim([0,1]);
%     set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
%         'TickDir','both','Box','on','FontName','Palatino Linotype',...
%         'FontSize',14);
%     ax=gca; obj=ax.findobj(); obj2=obj([13,10,7,4,12,9,6,3,11,8,5,2]);
%     legend(obj2,'Location','southoutside','Orientation','horizontal',...
%         'NumColumns',4,'FontSize',14,'Box','off','Interpreter','latex');
% end
% %
% cs={"L1","L2"}; 
% f=figure(3); clf; f.Color='w';  hold on
% f.Position=[450,380,675,500];
% for i=1:2
%     for j=1:nid
%         lnm="$\rm "+nms{j}+" - "+cs{i}+"$"; lwt=1;if j==1; lwt=2; end
%             p=plot(G*t{i,j},max(abs(dAv{i,j}),[],2),[clr{j},lst{i}],'DisplayName',lnm,...
%                 'LineWidth',lwt);
%     end
%     line([0,750],10.^[-6,-6],'Color','k','LineWidth',1,...
%         'LineStyle','--','HandleVisibility','off'); 
%     yscale('log'); ylim(10.^[-7 -1]);yticks(logspace(-7,-1,7)); xlim([0 750]);
%     xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',18); 
% 
%     % ylabel('$\rm max|\frac{\partial a_{ij}}{\partial t}|$',...
%     %     'Interpreter','latex','FontSize',14); 
%     ylabel("$log(err_{RK45})$",'Interpreter','latex','FontSize',18); 
%      grid on; 
%     set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
%         'TickDir','both','Box','on','FontName','Palatino Linotype',...
%         'FontSize',14);
%     legend('Box','off','Location','southoutside','Orientation','horizontal',...
%         'NumColumns',4,'Interpreter','latex');
% end
% %
% err=round(Avn1-Avn2,6)./round(Avn1,6)*100; err(isnan(err))=0;
% T=array2table(abs([err(:,[1 4 5],1) err(:,[1 4 5],2)]));
% T=varfun(@(x) num2str(x, '%.4f'),T);
% T.Properties.VariableNames=...
%     {'A_11_L1','A_22_L1','A_13_L1','A_11_L2','A_22_L2','A_13_L2'};
% T.Properties.RowNames={'RSC', 'FT', 'SRF','RPR'};
% %
% cs={"L1","L2"}; clst={'-', '-.'};
% f=figure(5); clf,f.Color='w'; hold on
% f.Position=[450,380,675,500];
% for i=1:2
%     for j=1:nid
%         lwt=1;if j==1; lwt=2; end
%         lnm="$\rm "+nms{j}+" - "+cs{i}+"$";
%         plot(Nerr{j,i},'Color',clr{j},'DisplayName',lnm,...
%             'Marker','x','MarkerSize',5,'LineStyle',clst{i},'LineWidth',lwt);
%     end
% end
% line([1,6],10.^[-6,-6],'Color','k','LineWidth',1.,...
%     'LineStyle','--','HandleVisibility','off'); 
% xlabel("$Iteration $",'FontSize',18,'Interpreter','latex');
% ylabel("$log(err_{NR})$",'FontSize',18,'Interpreter','latex');
% yscale('log');
% legend('Box','off','Location','southoutside','Orientation','horizontal',...
%         'NumColumns',4,'Interpreter','latex','FontSize',14);
% grid on; set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
%     'TickDir','both','Box','on','FontName','Palatino Linotype',...
%     'FontSize',14); xticks(1:6);yticks(logspace(-10,0,6));
% -------------------------------------------------------------------------
% Av0=[.3,1.e-10,1.e-10,.6,.1]'; Ar=1000; a=3;CLS=-4.1; G=1;
% dU=zeros(3); dU(6)=G;   
% CI={.0165;.0630;.0060}; ck={.988; .965; .900}; CM={.999; 1.010; .900}; 
% kpa={1/30; 1/30; 1/20}; ak=num2cell(1-[kpa{:}]'); bk={0.0;0.0;0.0};
% bta={[3.842 -17.86  525.  .1168  -5. ]*1e-4;
%      [37.28 -169.5  1750. -33.67 -100]*1e-4;
%      [4.643 -6.169  190.    9.65  7  ]*1e-4};
% for j=1:3, Dk{j,1}=diag([1, ck{j}, 1-ck{j}]); end
% val={[ak, bk, CI, CM], [ak, bk, CI, Dk], [kpa, bta]};
% id=[7.4, 7.3, 6.1] ;nid=length(id);
% Nerr=[];Avn2=[];Avn1=[]; Av=[]; t=[]; tn=[];
% for i=1:3
%     for j=1:nid
%         k=j+nid*(i-1);
%         var={Ar, a,@(t) dU,@(t) [],id(j), CLS,val{j}(i,:),{2,4}};
%         [Nerr{k}, Avn2(k,:)]=Newt([],Av0,var{:});[i,j]
%         [tn(k),Avn1(k,:),t{i,j},Av{i,j}]=Rung(T0,h,Tn,Av0,var{:});
%         for rk=1:length(t{i,j})
%             dAv{i,j}(rk,:)=ddA([],Av{i,j}(rk,:),v,var{2:end});
%         end
%     end
% end
% %
% Lns={'-', '--',':'}; pk=[1 4 5];
% nms={'iARD-RPR','pARD-RPR', 'iARD-RSC'}; clr={'r','g','c'};
% ttl={'40%wt.glass-fiber/PP', '31%wt.carbon-fiber/PP', '40%wt.glass-fiber/nylon'};
% for i=1:3
%     f=figure(i+3); clf; f.Color='w'; f.Position=[450,380,675,500];
%     grid on; hold on
%     % title(['\rm\it' ttl{i}]);
%     for j=1:nid
%         for k=1:3
%             Axx="$\rm "+nms{j}+"-a_{"+string(5*(pk(k)-k)+10+k)+"}$";
%             plot(G*t{i,j},Av{i,j}(:,pk(k)),Lns{k},'DisplayName',Axx,...
%                 'LineWidth',1.5,'Color',clr{j});
%         end
%     end
%     xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',18);
%     ylabel("$\rm a_{ij}$",'Interpreter','latex','FontSize',18);
%     set(gca,'TickDir','both','GridLineStyle','--','Box','on',...
%         'FontName','Palatino Linotype','FontSize',14);
%     xlim([0 G*Tn]); ylim([0,1]); 
%     legend('Location','southoutside','Orientation','horizontal','NumColumns',3,...
%         'FontSize',12,'Box','off','Interpreter','latex'); 
% end
% %
% f=figure(7); clf; f.Color='w'; f.Position=[450,380,675,500];
% grid on; hold on
% cs={"M1","M2","M3"}; clst={'-','--','-.'};
% for i=1:3
%     for j=1:nid
%         k=j+nid*(i-1); lnm="$\rm "+nms{j}+" - "+cs{i}+"$";
%         plot(G*t{i,j},max(abs(dAv{i,j}),[],2),Lns{i},'DisplayName',lnm,...
%             'LineWidth',1.5,'Color',clr{j});
%     end
% end
% line([0,3550],10.^[-6,-6],'Color','k','LineWidth',1,...
%     'LineStyle','--','HandleVisibility','off'); xlim([0,3550]);
%  yscale('log'); ylim(10.^[-7 -2]);yticks(logspace(-7,-2,6))
% xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',18); 
% % ylabel('$\rm max|\frac{\partial a_{ij}}{\partial t}|$',...
% %     'Interpreter','latex','FontSize',18); 
% ylabel("$log(err_{RK45})$",'Interpreter','latex','FontSize',14); 
% set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
%     'TickDir','both','Box','on','FontName',...
%     'Palatino Linotype','FontSize',14);
% legend('Location','southoutside','Orientation','horizontal',...
%     'NumColumns',3,'FontSize',11,'Box','off','Interpreter','latex');
% %
% err=round(Avn1-Avn2,6)./round(Avn1,6)*100; err(isnan(err))=0;
% T=array2table(abs(err(:,[1 4 5])));T=varfun(@(x) num2str(x, '%.4f'),T);
% T.Properties.VariableNames={'A_11','A_22','A_13'};
% Tn=['iARD-RPR';'pARD-RPR'; 'iARD-RSC'];nTn=size(Tn,1);
% Tn=[repmat(Tn,3,1),repelem(num2str((1:3)'),nTn,1)];nTn=size(Tn);
% T.Properties.RowNames=mat2cell(Tn,ones(nTn(1),1),nTn(2));
% %
% f=figure(8); clf; f.Color='w'; grid on; hold on
% f.Position=[450,380,675,500];
% for i=1:3
%     for j=1:nid
%         k=j+nid*(i-1);
%         lnm="$\rm "+nms{j}+" - "+cs{i}+"$";
%         plot(Nerr{k},clst{i},'Color',clr{j},'DisplayName',lnm,...
%             'Marker','x','MarkerSize',5,'LineWidth',1.5);
%     end
% end
% line([1,18],10.^[-6,-6],'Color','k','LineWidth',1.,...
%     'LineStyle','--','HandleVisibility','off'); yscale('log'); 
% xlabel("$Iteration $",'FontSize',18,'Interpreter','latex');
% ylabel("$log(err_{NR})$",'FontSize',18,'Interpreter','latex');
% grid on
% set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
%     'TickDir','both','Box','on','FontName',...
%     'Palatino Linotype','FontSize',14); 
% xlim([1,18]); ylim(10.^[-10 0]); yticks(logspace(-10,0,6))
% legend('Box','off','Location','southoutside','Orientation',...
%     'horizontal','NumColumns',3,'Interpreter','latex','FontSize',11);
%--------------------------------------------------------------------------
