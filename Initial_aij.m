clear,  close all, clc
[Rung, Newt, ddA]=Bfunc; warning('off');
Ar=1000; a=3; t0=0; ti=0.01; tn=500;
mkr={'o','+','*','.','x','_','|','s','d','^','v','>','<','p','h'};
clr=[
139, 69,19;0  ,255,  0;0   ,0  ,  0;255,192,205;0,255,  255;255,  0,255;
255,255, 0;255,165,   0;128,0  ,128;255,  0,  0;0,  0,  255]/255;
% -------------------------------------------------------------------------
v=(Ar^2-1)/(Ar^2+1); dU=zeros(3); dU(6)=1;
cls=-4.1; 
mdl_flg=1.0; mdl_par={1. .01}; % AT
% mdl_flg=2.0; mdl_par={.1 .01}; % RSC
% mdl_flg=3.3; c=.9868; D=[1,c,1-c] ; D=diag(D,0); mdl_par={.0169 D}; % pARD
% mdl_flg=3.2; mdl_par={.0562 .9977}; % iARD
% mdl_flg=7.0; mdl_par={0.9, .0, .01}; % AT-RPR
var={Ar, a, @(t) dU, @(t) [],mdl_flg,cls,mdl_par,{2,4}};
nj=11;
Aeps=1.e-10;A11=0.30; A23=0.10; A22=linspace(Aeps,1-A11-Aeps,nj);
for j=1:nj
    Av0=[A11,Aeps,Aeps,A22(j),A23]';
    [Nerr{j},Avn2(j,:)]=Newt([],Av0,var{:});
    [tn1(j) ,Avn1(j,:)]=Rung(t0,ti,tn,Av0,var{:});
end
err=round(Avn1-Avn2,6)./round(Avn1,6)*100  ;err(isnan(err))=0;
verr=vecnorm(err,2,2);
%
f=figure(4); clf,f.Color='w'; hold on
f.Position=[450,380,675,500];
for j=1:nj
    clrj=clr(j,:); lsty='-'; lwgt=1.5;
    if verr(j)>1, clrj=[128,128,128]/255; lsty='--'; lwgt=1.; end
    Axx=sprintf("$\\rm a_{22}^{0}=%.2f$",A22(j));
    plot(Nerr{j},'Color',clrj,'Marker','x','MarkerSize',5,...
        'LineStyle',lsty,'LineWidth',lwgt,'DisplayName',Axx);

end
line([1,8],10.^[-6,-6],'Color','k','LineWidth',1.,...
    'LineStyle','--','HandleVisibility','off'); 
xlabel("$Iteration$",'Interpreter','latex','FontSize',18); 
ylabel("$log(err_{NR})$",'Interpreter','latex','FontSize',18);
yscale('log'); legend('Box','off','Location','southwest',...
    'Interpreter','latex'); grid on
xlim([1,8]); ylim(10.^[-10,2]); yticks(logspace(-10,2,7));
set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
 'TickDir','both','Box','on','FontName','Palatino Linotype','FontSize',14);
%
lgd=legend; 
T=array2table(abs(err(:,[1 4 5],1)));T=varfun(@(x) num2str(x, '%.4f'),T);
T.Properties.RowNames=lgd.String;
T.Properties.VariableNames={'A_11','A_22','A_13'};