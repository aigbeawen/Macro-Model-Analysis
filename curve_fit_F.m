clear,  close, clc
[Rung, Newt, ddA]=Bfunc; Av0=[.4,1.e-4,1.e-4,.45,0.1]';
Ar=1000; a=3; t0=0; ti=.01; tn=100; v=(Ar^2-1)/(Ar^2+1);
mkr={'o','+','*','.','x','_','|','s','d','^','v','>','<','p','h'};
lgd={'$SS$','$SUA_1$','$SUA_2$','$UA$','$BA$','$PST_1$','$PST_2$',...
     '$SBA_1$','$SBA_2$'}; 
% -------------------------------------------------------------------------
id2=1.0; kp=1.; CI=0.01; val2={kp CI}; cls=-2.3; % -2.3, -3.1
% -------------------------------------------------------------------------
val1={1;[1,1;1,10]; 1; 1;[1,1;1,10];[1, 1;1, 10]};
nf=length(val1); jj=0;
for j=1:nf
    valj=val1{j}; ns(j)=size(valj,1);
    for k=1:ns(j)
        jj=jj+1; dU=DU(j,valj(k,:));
        var={Ar, a,@(t) dU , @(t) [],id2,cls,val2,{2,4}};
        [Nerr{jj}, Avn2(jj,:)                   ]=Newt([],Av0,var{:}); [j,k]
        [tnj(jj) , Avn1(jj,:),t{jj},Av{jj}]=Rung(t0,ti,tn,Av0,var{:}); tnj(jj)
        for rk=1:length(t{jj})
            dAv{jj}(rk,:)=ddA([],Av{jj}(rk,:),v,var{2:end});
        end
    end
end
lsty={'-','--','-.'}; clr={'r','g','g','b','c','b','m','m'};
% 
for m=1:3
    f=figure(m);clf;f.Color='w'; grid on; hold on
    n=-(m^2-6*m+4); ss=-16*m.^2+70*m-43; jj=0;
    Axx=['$\rm a_{' num2str(ss) '}$'];
    for j=1:nf
        nf2=ns(j);
        for k=1:nf2
            jj=jj+1;
            plot(t{jj},Av{jj}(:,n),'Color',clr{j},'LineStyle',lsty{k+(nf2>1)},...
                'LineWidth',.5,'DisplayName',lgd{jj});
        end
    end
    xlabel('\it\.{$\gamma$}t','Interpreter','latex','FontSize',18);
    ylabel(Axx,'Interpreter','latex','FontSize',18); xlim([0 70]); 
    legend('Location','best','Orientation','horizontal',...
        'Box','off', 'FontSize',14,'Interpreter','latex','NumColumns',3);
    ax=gca; obj=ax.findobj(); obj2=obj([1,10,5,4,7,9,8,6,3,2]);
    legend(obj2(2:end));   f.Position=[450,380,675,500];
    set(gca,'TickDir','both','GridLineStyle','--','Box','on',...
        'FontName','Palatino Linotype','FontSize',14);
end
%
err=round(Avn1-Avn2,5)./round(Avn1,5)*100        ;
err(isnan(err))=0; T=array2table(abs(err(:,[1 4 5]))) ;
T=varfun(@(x) num2str(x, '%.4f'),T);
T.Properties.VariableNames={'A_11','A_33','A_23'};
T.Properties.RowNames=lgd;
% -
f=figure(4);clf;f.Color='w'; grid on; hold on
jj=0;
for j=1:nf
    nf2=ns(j);
    for k=1:nf2
        jj=jj+1;
        plot(t{jj},max(abs(dAv{jj}),[],2),'Color',clr{j},'LineStyle',lsty{k+(nf2>1)},...
            'LineWidth',.5,'DisplayName',lgd{jj});
    end
end
line([0,70],10.^[-6,-6],'Color','k','LineWidth',1,...
    'LineStyle','--','HandleVisibility','off');
yscale('log'); ylim(10.^[-7 0]); yticks(logspace(-7,0,8)); xlim([0 70]); 
xlabel('\it\.{$\gamma$}t','Interpreter','latex');
% ylabel('$\rm max|\frac{\partial a_{ij}}{\partial t}|$',...
%     'Interpreter','latex','FontSize',14); 
ylabel("$log(err_{RK45})$",'Interpreter','latex','FontSize',18); 
legend('Location','best','Orientation','horizontal',...
    'Box','off', 'FontSize',14,'Interpreter','latex','NumColumns',3);
ax=gca; obj=ax.findobj(); obj2=obj([1,10,5,4,7,9,8,6,3,2]);
legend(obj2(2:end));   f.Position=[450,380,675,500];
set(gca,'TickDir','both','GridLineStyle','--','MinorGridLineStyle','none',...
    'Box','on', 'FontName','Palatino Linotype','FontSize',14);
%--------------------------------------------------------------------------
f=figure(5); clf; f.Color='w'; f.Position=[450,380,675,500];
grid on; hold on
jj=0;
for j=1:nf
    nf2=ns(j);
    for k=1:nf2
        jj=jj+1;
        plot(Nerr{jj},'Color',clr{j},'LineStyle',lsty{k+(nf2>1)},...
            'LineWidth',.5,'DisplayName',lgd{jj});
    end
end
xlabel("$Iteration $",'FontSize',18,'Interpreter','latex');
ylabel("$log(err_{NR})$",'FontSize',18,'Interpreter','latex');
line([1,8],10.^[-6,-6],'Color','k','LineWidth',1.,...
    'LineStyle','--','HandleVisibility','off'); 
yscale('log'); ylim(10.^[-10,2]); yticks(logspace(-10,2,7)); xlim([1,8]);
legend('Location','best','Orientation','horizontal',...
    'Box','off', 'FontSize',14,'Interpreter','latex','NumColumns',3);
ax=gca; obj=ax.findobj(); obj2=obj([1,10,5,4,7,9,8,6,3,2]);
    legend(obj2(2:end));    f.Position=[450,380,675,500];
set(gca,'GridLineStyle','--','MinorGridLineStyle','none',...
    'TickDir','both','Box','on',...
    'FontName','Palatino Linotype','FontSize',14);
% -------------------------------------------------------------------------
%
function dV=DU(stl,varargin)
    dat=varargin{:}; G=dat(1); GE=dat(end);
    E=G/GE; dV=zeros(3);I=eye(3);
    switch stl
        case 1 % Simple shear
            dV(6)=G;
        case 2 % Shearing/stretching
            dV([1,9,6,5])=[-E,-E, G, 2*E];
        case 3 % Uniaxial
            dV([1,9,5])=[-E,-E,2*E];
        case 4 % Biaxial
            dV([1,9,5])=[-2*E, E, E];
        case 5 % Shearing/planar stretching
            dV([1, 6, 5])=[-E,  G, E];
        case 6 % Balanced shear/biaxial elongation flow
            dV([1, 9, 6, 5])=[-2*E, E, G, E];
        case 7 % Triaxial
            dV([1,9,5])=[E, E, E];
        case 8 % Balanced shear/triaxial elongation flow
            dV([1,9,6,5])=[E, E, G, E];
    end
    dV=dV-trace(dV)/3*I; 
end
%%
function A=v2M(Av,flg)
for m=1:2
    for n=m:3
        k=2*(m-1)+n;
        switch flg
            case 1
                A(m,n)=Av(k); A(n,m)=A(m,n);
            case 2
                A(k,1)=Av(m,n);
        end
    end
end
if flg==1, A(3,3)=1-A(1,1)-A(2,2); end
end