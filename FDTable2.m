clear, clc
format shortg
[Rung, ~,ddA,dMdA]=BfuncC; v=1; a=3; re=1000; 
% A=rand(3); A=1/2*(A+A'); Ad=diag(A,0); Adev=.1*(A-diag(Ad,0));
% Aext=diag(Ad+1/3*(1-trace(A)),0);A=Aext+Adev;Av=[A(1,1:3), A(2,2:3)];
Av=[0.0622,0.0765,0.0398,0.5521,0.0186];
G=1; E=1; dU=zeros(3); dU([1,5, 6, 9])=[-2*E, E, G, E];
% -------------------------------------------------------------------------
mdl_id=[1.0,3.1:.1:3.3, 3.5,3.6, 5.0, 6.3, 7.2] ; 
mdl_nm={'FT','PT','iARD','pARD','WPT','Dz','NEM', 'pARD-RSC', 'iARD-RPR'};
cls_id=[1.1, 1.2, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 3.0,...
       -1.1,-2.1,-2.3,-2.4,-3.1,-3.2,-3.3,-3.4,-4.1,-4.2,-4.3];
cls_nm=["HYB_1","HYB_2","ISO","LIN","QDR","SF2","HL1","HL2",...
        "IBOF","ORS","ORT","NAT_1","ORW","NAT_2",...
        "WTZ","LAR32","ORW3","VST","FFLAR4","LAR4"];
n1=length(mdl_id); n2=length(cls_id); err=zeros(n1,n2);
%
delj=10.^(-12:1:-1); n3=length(delj); n4=6;
for j4=1:n4
    for j3=1:n3
        for j1=1:n1
            jmdl=mdl_id(j1);
            for j2=1:n2
                invars={v, a, @(t) dU ,@(t) [], jmdl,cls_id(j2),mdl(jmdl),{2,4}};
                [~, J1]=ddA([],Av, invars{:});
                [~ ,J2]=dMdA(@(x) ddA([], x, invars{:}), Av, delj(j3), j4);
                J2=reshape(J2,5,5); err(j1,j2)=norm(J1-J2);
            end
        end
        [emn(j3,j4),emx(j3,j4)]=bounds(err,'all'); eav(j3,j4)=mean(err,'all');
    end
end
fig1=figure(1); fig1.Color='w'; 
eb=errorbar(log10(delj),log10(eav(:,3)),...
    log10(eav(:,3)./emn(:,3)),log10(emx(:,3)./eav(:,3)));
eb.Marker='s'; eb.LineStyle='-.'; eb.LineWidth=1.5;
eb.MarkerEdgeColor='r'; eb.MarkerFaceColor='r'; eb.Color='k';
xlabel('$log(\delta) $','Interpreter','latex');
ylabel('$log(err_{2})$','Interpreter','latex');
grid('on');set(gca,'Tickdir','both','GridLineStyle','--','FontSize',12,...
    'FontName','Palatino Linotype');
%
fig2=figure(2); fig2.Color='w'; [eav2,kav]=min(eav,[],1); 
ind=(0:n4-1)'.*n3(:)+kav(:); emn2=emn(ind); emx2=emx(ind);
eb2=errorbar((1:n4)',log10(eav2(:)),...
    log10(eav2(:)./emn2(:)),log10(emx2(:)./eav2(:)));
eb2.Marker='s'; eb2.LineStyle='none'; eb2.LineWidth=1.5;
eb2.MarkerEdgeColor='r'; eb2.MarkerFaceColor='r'; eb2.Color='k';
labs={'BW_1','FW_1','CT_1','BW_2','FW_2','CT_2'};set(gca,'XTickLabel',labs);
ylabel('$log(err_{2})$','Interpreter','latex');
grid('on');set(gca,'Tickdir','both','GridLineStyle','--','FontSize',12,...
    'FontName','Palatino Linotype'); xlim([.95,6.05]);
% -------------------------------------------------------------------------
function f=mdl(flg)
flg1=floor(flg); flg2=round(10*(flg-flg1));
switch flg1
    case 1 % FT - SRF
        f={1. .0311}; % ['kp' 'CI']
    case 2 % RSC
        f={.2 .01}; % ['kp' 'CI']
end
switch flg1
    case 3
        switch flg2
            case 0 % IRD
                f={.0311}; % 'CI'
            case 1 % ARD
                f={[1.924,58.390,400.,.1168,0.]*1e-4}; % 'bta 1-5'
            case {2 4} % iARD
                f={.0562 .9977};     % ['CI' 'CM'] ;
            case 3 % pARD
                c=.9868; D=[1,c,1-c]; D=diag(D,0);
                f={.0169 D};      % ['CI' 'D' ] ;
            case 5 % WPT
                f={.0504 .9950}; % ['CI', 'w']
            case 6 % Dz
                n=[0,0,1];
                f={.0258, .0051, n}; % ['CI', 'Dz', 'n']
        end
    case 5 % NEM
        f={.01; .063};  % ['CI, U0']
    case 6 % pARD-RSC
        CI=.027; c=.95; kp=.8; D=[1,c,1-c]; D=diag(D,0);f={CI D};f=[{kp} f];
    case 7 % iARD-RPR
        CI=.063; CM=.995; lfa=.2; bta=.01;
        f={lfa bta, CI, CM}; 
end
end