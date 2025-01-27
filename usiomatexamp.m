addpath([pwd '/USIO/']) %mac
%%

%indnames=readtable('listUS.xlsx','ReadVariableNames',false);
indnames=readtable('IO_product_list.xlsx','Sheet', 'Sheet1','ReadVariableNames',false);
ind=~isnan(table2array(indnames(:,1)));
indnames=table2cell(indnames(ind,2));
%year = 2022;
year = 2020;

%IOT = readtable('REAL_USE.xlsx','Sheet',num2str(year),'ReadVariableNames',false);
%FU = readtable('REAL_FDAGG.xlsx','Sheet',num2str(year),'ReadVariableNames',false);
IOT = readtable('IOTUK.xlsx','Sheet',num2str(year),'ReadVariableNames',false);
FU = readtable('Aggreg_var.xlsx','Sheet',num2str(year),'ReadVariableNames',false);

%backlog = dlmread('Copy_of_iot_backlog.csv','\t',1,1);
%backlog = backlog(1:105); % retain 1:106

backlog = csvread('Copy_of_iot_backlog.csv',1,1)

IOT = table2array(IOT(2:end,2:end))'; indlist=1:size(IOT,2);
FU = table2array(FU(2:end,2:end));
b = FU(:,1); bg = sum(FU(:,[1 9 10 11]),2);
% total output 
pyi  = FU(:,20);
%final consumption
pfi11  = FU(:,1);
%labor income
empcomp =  FU(:,17) ;
% value added
va =  FU(:,19)  ;
% contruct the producer price index   cp: current price  pyp : previous year price

int_cons_CP = FU(:,21);
int_cons_pyp = FU(:,22);


disc=sum(IOT,2)==0; IOT(disc,:)=[]; IOT(:,disc)=[]; indlist(disc)=[]; b(disc)=[]; bg(disc)=[]; backlog(disc)=[]; pyi(disc)=[]; pfi11(disc)=[];empcomp(disc)=[];va(disc)=[];int_cons_CP(disc)=[];int_cons_pyp(disc)=[];
sig = IOT./sum(IOT,2); sig=sig(1:(end-1),1:(end-1)); indlist=indlist(1:end-1); n=size(sig,1);

% iteratively remove sectors that do not buy from or sell to other sectors.
% the while-loop ends after one iteration.
flag=1;
while flag
    ind=max([sum(sig,2)'==0;sum(sig,1)==0],[],1);sig(:,ind)=[];sig(ind,:)=[];b(ind)=[];bg(ind)=[];indlist(ind)=[];backlog(ind)=[];pyi(ind)=[];pfi11(ind)=[];empcomp(ind)=[];va(ind)=[];int_cons_CP(ind)=[];int_cons_pyp(ind)=[];
    if n==size(sig,1);flag=0;else n=size(sig,1);end
end
% note the set of sectors removed
misind=setdiff(1:size(FU,1),indlist); indnames(misind)

% we use private consumption as final demand; results are similar when including government spending.
b=b./sum(b); bg=bg./sum(bg);  I=eye(n); dw = b'*(eye(n)-sig)^-1;
corr(b,bg)
dlmwrite('indlist.csv',[indlist']);

% calculate Antras et al. upstreamness measure
th=sig.*repmat(dw',1,n)./repmat(dw,n,1); up=sum((eye(n)-th)^-1,1);

rho = 0.04; del=3.2/12;
rhodel = rho*del; r=1/(1+rhodel);
%[evr,ev]=eigs(sig,n,'largestabs');
[evr,ev]=eig(sig);
[~,ind]=sort(abs(diag(ev)),'descend');
ev=ev(ind,ind);
evr=evr(:,ind);
% order the first 4 eigenvectors so that the majority entries are positive
for i=2:4
    if mean(evr(:,i))<0
        evr(:,i)=-evr(:,i);
    end
end
evl=evr^-1;
lamb = diag(ev); 
ev_linv_dw=diag(1./(1-lamb));
ev_linv_v=diag(lamb./(1-lamb).*(1-r)./(1-r*lamb));

sigrm=sig./(1+rho*repmat(backlog',n,1)/12);
v=b'*((I-sig)^-1 - (I-sig./(1+rhodel))^-1);
vH=b'*((I-sig)^-1 - (I-sigrm)^-1);
eta=v./dw;
etaH=vH./dw;

[~,ind]=sort(up,'descend'); th=th(ind,ind);
dlmwrite('usa_th.csv',[indlist',th]);
inf=real((diag(b'*evr)*ev_linv_v*evl)'); cuminf=cumsum(inf,2);
dlmwrite('eigenreal_usa_cuminf.csv',[indlist',cuminf(:,[1:6,end])]);
dwinf = real((diag(b'*evr)*ev_linv_dw*evl)'); cumdwinf=cumsum(dwinf,2);

disp('Table 4: v')
R2=[];indvec=[];coefvec=[];
for i=[1:5 10 15]
    [coef,~,~,~,stats]=regress(cuminf(:,end),[ones(n,1), cuminf(:,i)]); %[i coef(2),stats(1)]
    indvec=[indvec,i];
    coefvec=[coefvec,coef(2)];
    R2=[R2,stats(1)];
end
[indvec;coefvec;R2]

disp('Table 5: domar weight')
R2=[];indvec=[];coefvec=[];
for i=[1:5 10 15]
%for i=[166:172]
    [coef,~,~,~,stats]=regress(cumdwinf(:,end),[ones(n,1), cumdwinf(:,i)]); %[i coef(2),stats(1)]
    indvec=[indvec,i];
    coefvec=[coefvec,coef(2)];
    R2=[R2,stats(1)];
end
[indvec;coefvec;R2]

%%
disp('Table 1')
[~,ind]=sort(vH,'descend');
indnames(indlist(ind(1:10)))
indnames(indlist(ind(end:-1:(end-9))))

disp('Table 2')
[~,ind]=sort(dw,'descend');
indnames(indlist(ind(1:10)))
[~,ind]=sort(etaH,'descend');
indnames(indlist(ind(1:10)))

disp('Table 3: Pearson and Spearman correlations')
corr([etaH',eta',up'])
corr([etaH',eta',up'],'Type','Spearman')

%%
disp('Table OA.1')
[~,ind]=sort(abs(evr(:,1)),'descend');
names=indnames(indlist(ind(1:10)));
[evr_display,ind]=sort(evr(ind(1:10),1),'descend');
names(ind)

[~,ind]=sort(abs(evr(:,2)),'descend');
names=indnames(indlist(ind(1:10)));
[evr_display,ind]=sort(evr(ind(1:10),2),'descend');
names(ind)

[~,ind]=sort(abs(evr(:,3)),'descend');
names=indnames(indlist(ind(1:10)));
[evr_display,ind]=sort(evr(ind(1:10),3),'descend');
names(ind)


%% figure 5
Q=@(t) -b'*sig*(I-sig)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-sig)*t);
tvec=0:.01:5;
qvec=zeros(n,length(tvec));
for i=1:length(tvec)
     qvec(:,i)=(Q(tvec(i)));
end
%qvec=-qvec./repmat(qvec(:,1),1,length(tvec)); 

qvec=qvec*100; % this was multiplied by 100 in LT paer
halflives=zeros(n,1);
for i=1:n
    halflives(i)=find(qvec(i,:)>-50,1);
end
% one unit of time is a year; multiply the index by 12/100 to convert to months
halflives=halflives*0.12;
%%
try close(fig1); catch; end
ps=100;
fig1=figure(1);
ax=axes; 
hold on; 
xlim([0 max(tvec)]);
ylim([-5 0.1]);

%   {'Construction'                                                   }
%    {'Financial services, except insurance and pension funding'       }
%    {'Computer programming, consultancy and related services'         }
%    {'Employment services'                                            }
%    {'Electricity, transmission and distribution'                     }
%    {'Services auxiliary to financial services and insurance services'}
%    {'Computer, electronic and optical products'                      }
%    {'Real estate services on a fee or contract basis'                }
%    {'Services of head offices; management consulting services'       }
%    {'Advertising and market research services'                       }

%legs = {"Warehousing and support services for transportation" ,"Residential Care and Social Work Activities",...
%"Dyestuffs, agro chemicals", "Accounting, bookkeeping and auditing services; tax consulting services",...
%"Advertising and market research services",...
%"Services auxiliary to financial services and insurance services"}

legs = {"Construction" ,"Financial services, except insurance and pension funding",...
"Computer programming, consultancy and related services",...
"Employment services",...
"Electricity, transmission and distribution",...
"Services auxiliary to financial services and insurance services"}

%
kvec=[57 72 70 86 51 74];

%kvec=[9 10 6 7 8 5];

%line([0,0],[-1,0],'linewidth',1,'linestyle','--','color','black','HandleVisibility','off')
bl=0.2*ones(1,3); gr=0.7*ones(1,3);
plot(tvec,qvec(kvec(1),:),'linewidth',2,'linestyle','-','color','#337316');
plot(tvec,qvec(kvec(2),:),'linewidth',2,'linestyle',':','color','#055BA6');
plot(tvec,qvec(kvec(3),:),'linewidth',2,'linestyle','--','color','#023E73');
plot(tvec,qvec(kvec(4),:),'linewidth',2,'linestyle','-','color','#223F59');
plot(tvec,qvec(kvec(5),:),'linewidth',2,'linestyle',':','color','#F25922');
plot(tvec,qvec(kvec(6),:),'linewidth',2,'linestyle','--','color','#D92344');
%plot([-tvec((ps+1):0:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');
plot([-tvec((ps+1):-1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');

h=legend(legs);
set(h,'FontSize',8,'Location','east');
set(ax,'FontSize', 16, 'xtick',[0 1 2 3 4 5]);
set(ax,'FontSize', 16, 'ytick',[-5 -4 -3 -2 -1 0]);
text(max(tvec)-10,-1.15,"time",'FontSize',18)
xl=xlabel("Years since shock",'FontSize',18,'fontweight','normal');
ylabel(["GDP loss", "(relative to time-zero impact)"],'FontSize',18,'fontweight','normal');
set(gcf,'OuterPosition', [800, 800, 700, 550]);
exportgraphics(fig1,'gdploss_by_sector.pdf');

%%
legs = {"Computer, electronic and optical products","Real estate services on a fee or contract basis"...
"Services of head offices; management consulting services"...
"Advertising and market research services"}

kvec=[39 75 79 82];
try close(fig1); catch; end
ps=100;
fig1=figure(1);ax=axes; hold on; 
xlim([0 max(tvec)]); hold on; 
ylim([-4.5 0.1]);  hold on;
%line([0,0],[-120,0],'linewidth',1,'linestyle','--','color','black','HandleVisibility','off')
bl=0.2*ones(1,3); 
gr=0.7*ones(1,3);
plot(tvec,qvec(kvec(1),:),'linewidth',2,'linestyle','-','color','#055BA6');
plot(tvec,qvec(kvec(2),:),'linewidth',2,'linestyle',':','color','#023E73');
plot(tvec,qvec(kvec(3),:),'linewidth',2,'linestyle','--', 'color','#223F59');
plot(tvec,qvec(kvec(4),:),'linewidth',2,'linestyle','-.','color','#D92344');
plot([-tvec((ps+1):-1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');


%plot([tvec((ps+1):1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');
h=legend(legs);
set(h,'FontSize',8,'Location','best');
set(ax,'FontSize', 16, 'xtick',[0 1 2 3 4 5]);
set(ax,'FontSize', 16, 'ytick',[-4 -3 -2 -1 0]);

%ytickformat('percentage');
text(max(tvec)-10,-1.15,"time",'FontSize',18)
xl=xlabel("Years since shock",'FontSize',18,'fontweight','normal');
ylabel(["GDP loss", "(relative to time-zero impact)"],'FontSize',18,'fontweight','normal');
set(gcf,'OuterPosition', [800, 800, 700, 550]);
exportgraphics(fig1,'gdploss_by_sector2.pdf');


%%

kvec=[57 72 70 86 51 74 39 75 79 82];

disp(['half lives'])
mean(halflives)
mean(halflives(kvec(1:3)))
mean(halflives(kvec(4:6)))

disp(['impact after one year'])
mean(qvec(kvec(1:3),100))
mean(qvec(kvec(4:6),100))

disp(['relative impact after three years'])
mean(qvec(kvec(4:6),300))/mean(qvec(kvec(1:3),300))

%%  Barrot model extension 
Sig_E = sig;
pfi  = pfi11';

%%  IO matrix
pxij = sig;
pmXi = sum(pxij,2)

%%
% number of sectors
N = length(pxij);
% intermediate good matrix (coefficient)
Omega = diag( pmXi.^(-1) ) * pxij;







%%
%%% eta_i = value added share (of revenue) in sector i
%eta = 1 - pmXi./pyi' ;
eta = 1 - pmXi./pyi ;
    
%%% phi_i = share of final demand i in total output of i
phi = pfi ./ pyi;

%%
%%% Delta_ij =  expenditure on j by i as a share of total production of j
Delta = pxij * diag( pyi.^(-1) ) ;
    
%%% psi_i = share of final demand i in to total final demand
psi = pfi ./sum( pfi);

%%  matrix share of expenduture on j by i in total production j
%normalized the matrix
Gamma=zeros(size(pxij));
for i=1:N
Gamma(i,:) = pxij(i,:) ./ pyi(i);
end

%%  
    
%Labor income share in Value Added 
    
wli_over_vhi = (empcomp./va)';       



%%

    
%% steps to compute the Intermediate good price
% σ is the (constant) elasticity 
% of substitution across goods.
% θ is the (constant) elasticity of substitution between the capital/labor 
% bundle and the intermediate input bundle, 
% and ε is the (constant) elasticity of substitution 
% between intermediate inputs
sigma =  0.9; 
theta =  0.5;
epsilon= 0.001;  
%%

% # Check for zero values before division

int_cons_pyp_nonzero = int_cons_pyp;
int_cons_pyp_nonzer(int_cons_pyp == 0) = 1;% Replace zero with

int_cons_CP_nonzero = int_cons_CP;
int_cons_CP_nonzero(int_cons_CP == 0)= 1 ; % Replace zero with
%

inf_V = int_cons_CP_nonzero./int_cons_pyp_nonzero;

%  assume that:  price_{t-1} was one then price_{t} = inflation_{t}
%



% we can use either the price index of producers
% power of prices
pi_hat = inf_V ;

%pi_hat = [2 3 2 2 3 4 6 2 2 2 3]; % we can use either the price index of producers
%power of prices
pi_hat_1MoinsEpsi = pi_hat.^(1-epsilon);
pi_hat_1MoinsThetai  = pi_hat.^(1-theta);
    
    
%Intermediate Bundle Price
%        Pmi_hat = exp( Omega * log(pi_hat)  );
Pmi_hat = ( Omega .* pi_hat_1MoinsEpsi ).^(1/(1-epsilon));
%%

%%% gamma_i = labor income share in value added in sector i
gamma = wli_over_vhi;

% labor and capital 
li_ai = 0.3 + (1-0.3).*rand(n,1);
ki_ai =  ones(n,1);

%Value added bundle quantity 
ai = (li_ai').^gamma .* (ki_ai').^(1-gamma);
% value added share in total output 
%neta = va/pyi;

%wli_over_GDP
lambda = wli_over_vhi.*va./sum(pfi);
%rli_over_GDP
rho = (1-wli_over_vhi).*va./sum(pfi);
vv = va./ai;


%%

% Matrix I 
I=eye(n); 


% The parameter d captures the ease of adjustment when input expands
d=0.12; % 12

% data  b is the consumption vector Final consumption expenditure) 
%b = [30 120 90 80 160 80 160 160 300 300 310]'; 
%b = b;  pfi;

b_vec = b; %./sum(b);


shock = +1; % negative or positive
% Update the specification of output
% Q is output (industries intermediate consumption  - intermediate inputs mij ) function 
%Q1=@(t)  ( b_vec' * shock*Sig_E*(I-Sig_E)^-1*expm(-1/d*(I-Sig_E)*t)  )  ;
Q1=@(t)  ( b_vec' * shock*Sig_E*(I-Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t)  )  ;

%% Pmi_hat 
Q2=@(t)  real( ( ( pi_hat_1MoinsEpsi' * (shock *Sig_E*(I-Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t) ) )  ).^(1/(1-epsilon)) );
%% Pi_hat
%Q4=@(t)  real( (  shock^(theta-1) * ( eta'.*ai'.^(1-theta) + (1-eta)'.*( ( pi_hat_1MoinsEpsi'  * (Sig_E*(I-Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t)  )  ).^(1/(1-epsilon))).^(1-theta) ) ).^(1/(1-theta)) );
Q4=@(t)  real( (  shock^(theta-1) * (  eta'.*ai.^(1-theta) +  (1-eta)'.*( ( pi_hat_1MoinsEpsi'  * (Sig_E*(I-Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t)  )  ).^(1/(1-epsilon))).^(1-theta) ) ).^(1/(1-theta)) )
%% pi_hat alterantive (price of final good)
Q5=@(t)  real( (  eta'*(vv).^((theta-1)/theta) + (1-eta)'*( Pmi_hat *shock*Sig_E*(I-Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t)  ).^((theta-1)/theta)  ).^(theta/(theta-1)) );
%% alternative Pi_hatx
P = 1;
Pi_hat_Alt= ( P.^(1-sigma)./(psi) ).^(1./(1-sigma));

%Q3=@(t)  (  eta*(ai)^((theta-1)/theta) + (1-eta)*( Pmi_hat *shock*Sig_E*(I-Sig_E)^-1*expm(-1/d*(I-Sig_E)*t)  )^((theta-1)/theta)  )^(theta/(theta-1)) ;
%%

% time dimension 20 years
%tvec=0:1:20;  
tvec=0:.01:5;

% vector of sectorl output (3 sectors)
qvec1=zeros(n,length(tvec));


for i=1:length(tvec)
     qvec1(:,i)=(Q1(tvec(i)));
end
%%

% time dimension 20 years
%tvec=0:1:20;  

% vector of sectorl output (3 sectors)
qvec2=zeros(n,length(tvec));


for i=1:length(tvec)
     qvec2(:,i)=(Q2(tvec(i)));
end

%%
% vector of sectorl output (3 sectors)
qvec4=zeros(n,length(tvec));


for i=1:length(tvec)
     qvec4(:,i)=(Q4(tvec(i)));
end

%%
% vector of sectorl output (3 sectors)
qvec5=zeros(n,length(tvec));


for i=1:length(tvec)
     qvec5(:,i)=(Q5(tvec(i)));
end

%%

%qvec1=qvec1./repmat(qvec1(:,1),1,length(tvec)); 
qvec1=qvec1*100; % this was multiplied by 100 in LT paer
%%

try close(fig1); catch; end
ps=100;
fig11=figure(1);ax=axes; hold on; xlim([0 max(tvec)]);ylim([-0.1 5]);

legs = {"Construction" ,"Financial services, except insurance and pension funding",...
"Computer programming, consultancy and related services",...
"Employment services",...
"Electricity, transmission and distribution",...
"Services auxiliary to financial services and insurance services"}

%
kvec=[57 72 70 86 51 74];

%line([0,0],[-1,0],'linewidth',1,'linestyle','--','color','black','HandleVisibility','off')
bl=0.2*ones(1,3); gr=0.7*ones(1,3);
plot(tvec,qvec1(kvec(1),:),'linewidth',2,'linestyle','-','color','#337316');
plot(tvec,qvec1(kvec(2),:),'linewidth',2,'linestyle',':','color','#055BA6');
plot(tvec,qvec1(kvec(3),:),'linewidth',2,'linestyle','--','color','#023E73');
plot(tvec,qvec1(kvec(4),:),'linewidth',2,'linestyle','-','color','#223F59');
plot(tvec,qvec1(kvec(5),:),'linewidth',2,'linestyle',':','color','#F25922');
plot(tvec,qvec1(kvec(6),:),'linewidth',2,'linestyle','--','color','#D92344');
%plot([-tvec((ps+1):-1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');
plot([-tvec((ps+1):-1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');

h=legend(legs);
set(h,'FontSize',8,'Location','east');
set(ax,'FontSize', 16, 'xtick',[0 1 2 3 4 5], 'ytick',[0 1  2  3 4 5]);
%ytickformat('percentage');
text(max(tvec)-10,-1.15,"time",'FontSize',18)
xl=xlabel("Years since shock",'FontSize',18,'fontweight','normal');
ylabel(["GDP Gain", "(relative to time-zero impact)"],'FontSize',18,'fontweight','normal');
set(gcf,'OuterPosition', [800, 800, 700, 550]);
exportgraphics(fig11,'gdpgain_by_sector2.pdf');




%%
%try close(fig11); catch; end
%ps=100;
fig2=figure(2);
ax=axes; hold on; 
xlim([0 max(tvec)]);
ylim([-0.1 5]);


legs = {"Construction" ,"Financial services, except insurance and pension funding",...
"Computer programming, consultancy and related services",...
"Employment services",...
"Electricity, transmission and distribution",...
"Services auxiliary to financial services and insurance services"}

%
kvec=[57 72 70 86 51 74];


plot(tvec,qvec2(kvec(1),:),'linewidth',2,'linestyle','-','color','#337316');
plot(tvec,qvec2(kvec(2),:),'linewidth',2,'linestyle',':','color','#055BA6');
plot(tvec,qvec2(kvec(3),:),'linewidth',2,'linestyle','--','color','#023E73');
plot(tvec,qvec2(kvec(4),:),'linewidth',2,'linestyle','-','color','#223F59');
plot(tvec,qvec2(kvec(5),:),'linewidth',2,'linestyle',':','color','#F25922');
plot(tvec,qvec2(kvec(6),:),'linewidth',2,'linestyle','--','color','#D92344');
plot([-tvec((ps+1):-1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');
h=legend(legs);
set(h,'FontSize',8,'Location','east');
set(ax,'FontSize', 16, 'xtick',[0 1 2 3 4 5], 'ytick',[0 1 2 3 4 5]);
%ytickformat('percentage');
text(max(tvec)-10,-1.15,"time",'FontSize',18)
xl=xlabel("Years since shock",'FontSize',18,'fontweight','normal');
ylabel(["Price of Intermediate good", "(relative to time-zero impact)"],'FontSize',18,'fontweight','normal');
set(gcf,'OuterPosition', [800, 800, 700, 550]);

exportgraphics(fig2,'priceintermediatgood_by_sector.pdf');

%%

hold off;
%%
fig8=figure(10);
ax=axes; hold on; 
xlim([0 max(tvec)]);
ylim([0 3]);

legs = {"Construction" ,"Financial services, except insurance and pension funding",...
"Computer programming, consultancy and related services",...
"Employment services",...
"Electricity, transmission and distribution",...
"Services auxiliary to financial services and insurance services"}

%
kvec=[57 72 70 86 51 74];

plot(tvec,qvec5(kvec(1),:),'linewidth',2,'linestyle','-','color','#337316');
plot(tvec,qvec5(kvec(2),:),'linewidth',2,'linestyle',':','color','#055BA6');
plot(tvec,qvec5(kvec(3),:),'linewidth',2,'linestyle','--','color','#023E73');
plot(tvec,qvec5(kvec(4),:),'linewidth',2,'linestyle','-','color','#223F59');
plot(tvec,qvec5(kvec(5),:),'linewidth',2,'linestyle',':','color','#F25922');
plot(tvec,qvec5(kvec(6),:),'linewidth',2,'linestyle','--','color','#D92344');
h=legend(legs);
set(h,'FontSize',8,'Location','east');
set(ax,'FontSize', 16, 'xtick',[0 1 2 3 4 5], 'ytick',[0 1 2 3 ]);
%ytickformat('percentage');
text(max(tvec)-10,-1.15,"time",'FontSize',18)
xl=xlabel("Years since shock",'FontSize',18,'fontweight','normal');
ylabel(["Price of Final good", "(relative to time-zero impact)"],'FontSize',18,'fontweight','normal');
set(gcf,'OuterPosition', [800, 800, 700, 550]);

exportgraphics(fig8, 'figalt1.pdf'); % Save as PNG format

hold off

%%  vH dw





fig51=figure(51);


% Generate x values for the line  impact
x_values = vH *100

% Calculate corresponding y values using the fitted coefficients
y_values = dw

% Plot the linear regression line
%plot(x_values, y_values, 'LineWidth', 3, 'Color', '#FFA500'); % Plot the linear regression line
% Fit a linear regression model
mdl = fitlm(x_values, y_values);

% Generate values for the regression line
x_reg = linspace(min(x_values), max(x_values), 100);
y_reg = predict(mdl, x_reg');

% Plot the scatter plot
scatter(x_values, y_values,'Color','#0F95D7'); hold on
%scatter(kappae_v_var, output_var,'Color', '#0F95D7'); hold on
% Plot the linear regression line
plot(x_reg, y_reg, 'LineWidth', 1.5, 'Color', '#D92344'); % Plot the linear regression line

% Add zero lines and text annotations
xline(0, '--k', 'LineWidth', 1); % Vertical line at x = 0
yline(0, '--k', 'LineWidth', 1); % Horizontal line at y = 0

ax = gca;
ax.FontSize = 18;
ax.XGrid = 'off';
ax.YGrid = 'on';
axis tight; 
xlabel('Impact','Interpreter','latex'); % Add labels and title
ylabel('Domar weight','Interpreter','latex');


ylim([-0.003 0.1])
xlim([-0.003 0.1])

exportgraphics(fig51, 'fig6.pdf'); % Save as PNG format
hold off

%%   figure 6

% disp('Table 6')
%[~,ind]=sort(dw,'descend');
%indnames(indlist(ind(1:10)))
%[~,ind]=sort(etaH,'descend');
%indnames(indlist(ind(1:10)))


Sorted_sig = sig;


%%


% Define the RGB values for the two colors
color2 = [217, 4, 4] / 255; 
color1 = [242, 242, 242] / 255;

% Create a colormap interpolating between the two colors
custom_colormap = [linspace(color1(1), color2(1), 100)', ...
                   linspace(color1(2), color2(2), 100)', ...
                   linspace(color1(3), color2(3), 100)'];


fig5=figure(5);

imagesc(Sorted_sig);
colormap(custom_colormap);
h= colorbar; 
h.Ticks = [];


axis equal; 
axis tight; 


%saveas(gcf, 'figLiu6.pdf'); % Save as PNG format
exportgraphics(fig5,'figLiu6.pdf');
%%

color1 = [242, 242, 242] / 255; % Light gray
color2 = [1, 28, 63]/255; % Red

% Number of colors in the colormap
nColors = 256;

% Create the custom colormap
customMap = [linspace(color1(1), color2(1), nColors)', ...
             linspace(color1(2), color2(2), nColors)', ...
             linspace(color1(3), color2(3), nColors)'];

fig6=figure(6);
customYTicks = [1, 20, 40, 60, 80, 100]; % Example tick positions, adjust as needed
customYTickLabels = fliplr(customYTicks); % Labels in descending order

transposedMatrix = Sorted_sig;
flippedTransposedMatrix = flipud(transposedMatrix);



contour(flippedTransposedMatrix); % or contour(Z) if X and Y are not needed
colorbar; % Adding a colorbar
colormap(customMap);
%title('Contour Plot of Matrix Z');
xlabel(['Suppliers']);
ylabel('Buyers');
yticks(customYTicks); % Set custom y-axis ticks
yticklabels(customYTickLabels); % Set custom y-axis labels

exportgraphics(fig6,'figcountourLiu6.pdf');

%%

dlmwrite('figure_6.csv',[Sorted_sig]);

%%

%%  targeted policy  (manufaturing sector) referenvc year 2020 
% subsidy for firms
%Sub= 1.15;

Sub_vec = ones(1, 101);

% Step 2: Set the specified ranges to zeros
vec(18:20) = 0.5;  % Set rows 18 to 24 to zero
vec(23:50) = 0.5;  % Set rows 30 to 35 to zero

% Display the result
disp(vec);

%18	Textiles               	Manufacturing	2.367491	CPA_C13
%19	Wearing apparel              	Manufacturing	2.367491	CPA_C14
%20	Leather and related products            	Manufacturing	2.367491	CPA_C15
%23	Printing and recording services            	Manufacturing	2.367491	CPA_C18
%24	Coke and refined petroleum products           	Manufacturing	2.367491	CPA_C19
%25	Paints, varnishes and similar coatings, printing ink and mastics       	Manufacturing	2.367491	CPA_C203
%26	Soap and detergents, cleaning and polishing preparations, perfumes and toilet preparations     	Manufacturing	2.367491	CPA_C204
%27	Other chemical products             	Manufacturing	2.367491	CPA_C205
%28	Industrial gases, inorganics and fertilisers (all inorganic chemicals) - 20.11/13/15      	Manufacturing	2.367491	CPA_C20A
%29	Petrochemicals - 20.14/16/17/60             	Manufacturing	2.367491	CPA_C20B
%30	Dyestuffs, agro-chemicals - 20.12/20            	Manufacturing	2.367491	CPA_C20C
%31	Basic pharmaceutical products and pharmaceutical preparations          	Manufacturing	2.367491	CPA_C21
%32	Rubber and plastic products            	Manufacturing	2.367491	CPA_C22
%33	Cement, lime, plaster and articles of concrete, cement and plaster 	Manufacturing	2.367491	CPA_C235_6
%34	Glass, refractory, clay, other porcelain and ceramic, stone and abrasive products - 23.1-4/7-9   	Manufacturing	2.367491	CPA_C23OTHER
%35	Basic iron and steel            	Manufacturing	2.367491	CPA_C241_3
%36	Other basic metals and casting           	Manufacturing	2.367491	CPA_C244_5
%37	Weapons and ammunition             	Manufacturing	2.367491	CPA_C254
%38	Fabricated metal products, excl. machinery and equipment and weapons & ammunition - 25.1-3/25.5-9   	Manufacturing	2.367491	CPA_C25OTHER
%39	Computer, electronic and optical products           	Manufacturing	2.367491	CPA_C26
%40	Electrical equipment              	Manufacturing	2.367491	CPA_C27
%41	Machinery and equipment n.e.c.            	Manufacturing	2.367491	CPA_C28
%42	Motor vehicles, trailers and semi-trailers           	Manufacturing	2.367491	CPA_C29
%43	Ships and boats             	Manufacturing	2.367491	CPA_C301
%44	Air and spacecraft and related machinery          	Manufacturing	2.367491	CPA_C303
%45	Other transport equipment - 30.2/4/9           	Manufacturing	2.367491	CPA_C30OTHER
%46	Furniture               	Manufacturing	2.367491	CPA_C31
%47	Other manufactured goods             	Manufacturing	2.367491	CPA_C32
%48	Repair and maintenance of ships and boats         	Manufacturing	2.367491	CPA_C3315
%49	Repair and maintenance of aircraft and spacecraft         	Manufacturing	2.367491	CPA_C3316
%50	Rest of repair; Installation - 33.11-14/17/19/20          	Manufacturing	2.367491	CPA_C33OTHER


%%

indnames=readtable('IO_product_list.xlsx','Sheet', 'Sheet1','ReadVariableNames',false);
ind=~isnan(table2array(indnames(:,1)));
indnames=table2cell(indnames(ind,2));
year = 2020;


IOT = readtable('IOTUK.xlsx','Sheet',num2str(year),'ReadVariableNames',false);
FU = readtable('Aggreg_var.xlsx','Sheet',num2str(year),'ReadVariableNames',false);

%backlog = dlmread('Copy_of_iot_backlog.csv','\t',1,1);
%backlog = backlog(1:105); % retain 1:106

backlog = csvread('Copy_of_iot_backlog.csv',1,1)

IOT = table2array(IOT(2:end,2:end))'; indlist=1:size(IOT,2);
FU = table2array(FU(2:end,2:end));
b = FU(:,1) ; bg = sum(FU(:,[1 9 10 11]),2);
% total output 
pyi  = FU(:,20) ;
%final consumption
pfi11  = FU(:,1)  ;
%labor income
empcomp =  FU(:,17) ;
% value added
va =  FU(:,19);

% contruct the producer price index   cp: current price  pyp : previous year price

int_cons_CP = FU(:,21);
int_cons_pyp = FU(:,22);


disc=sum(IOT,2)==0; IOT(disc,:)=[]; IOT(:,disc)=[]; indlist(disc)=[]; b(disc)=[]; bg(disc)=[]; backlog(disc)=[]; pyi(disc)=[]; pfi11(disc)=[];empcomp(disc)=[];va(disc)=[];int_cons_CP(disc)=[];int_cons_pyp(disc)=[];
sig = IOT./sum(IOT,2); sig=sig(1:(end-1),1:(end-1)); indlist=indlist(1:end-1); n=size(sig,1);

% iteratively remove sectors that do not buy from or sell to other sectors.
% the while-loop ends after one iteration.
flag=1;
while flag
    ind=max([sum(sig,2)'==0;sum(sig,1)==0],[],1);sig(:,ind)=[];sig(ind,:)=[];b(ind)=[];bg(ind)=[];indlist(ind)=[];backlog(ind)=[];pyi(ind)=[];pfi11(ind)=[];empcomp(ind)=[];va(ind)=[];int_cons_CP(ind)=[];int_cons_pyp(ind)=[];
    if n==size(sig,1);flag=0;else n=size(sig,1);end
end
% note the set of sectors removed
misind=setdiff(1:size(FU,1),indlist); indnames(misind)

% we use private consumption as final demand; results are similar when including government spending.
b=b./sum(b); bg=bg./sum(bg);  I=eye(n); dw = b'*(eye(n)-sig)^-1;
corr(b,bg)
dlmwrite('indlist.csv',[indlist']);

% calculate Antras et al. upstreamness measure
th=sig.*repmat(dw',1,n)./repmat(dw,n,1); up=sum((eye(n)-th)^-1,1);

rho = 0.04; del=3.2/12;
rhodel = rho*del; r=1/(1+rhodel);
%[evr,ev]=eigs(sig,n,'largestabs');
[evr,ev]=eig(sig);
[~,ind]=sort(abs(diag(ev)),'descend');
ev=ev(ind,ind);
evr=evr(:,ind);
% order the first 4 eigenvectors so that the majority entries are positive
for i=2:4
    if mean(evr(:,i))<0
        evr(:,i)=-evr(:,i);
    end
end
evl=evr^-1;
lamb = diag(ev); 
ev_linv_dw=diag(1./(1-lamb));
ev_linv_v=diag(lamb./(1-lamb).*(1-r)./(1-r*lamb));

sigrm=sig./(1+rho*repmat(backlog',n,1)/12);
v=b'*((I-sig)^-1 - (I-sig./(1+rhodel))^-1);
vH=b'*((I-sig)^-1 - (I-sigrm)^-1);
eta=v./dw;
etaH=vH./dw;

[~,ind]=sort(up,'descend'); th=th(ind,ind);
dlmwrite('usa_th.csv',[indlist',th]);
inf=real((diag(b'*evr)*ev_linv_v*evl)'); cuminf=cumsum(inf,2);
dlmwrite('eigenreal_usa_cuminf.csv',[indlist',cuminf(:,[1:6,end])]);
dwinf = real((diag(b'*evr)*ev_linv_dw*evl)'); cumdwinf=cumsum(dwinf,2);

disp('Table 4: v')
R2=[];indvec=[];coefvec=[];
for i=[1:5 10 15]
    [coef,~,~,~,stats]=regress(cuminf(:,end),[ones(n,1), cuminf(:,i)]); %[i coef(2),stats(1)]
    indvec=[indvec,i];
    coefvec=[coefvec,coef(2)];
    R2=[R2,stats(1)];
end
[indvec;coefvec;R2]

disp('Table 5: domar weight')
R2=[];indvec=[];coefvec=[];
for i=[1:5 10 15]
%for i=[166:172]
    [coef,~,~,~,stats]=regress(cumdwinf(:,end),[ones(n,1), cumdwinf(:,i)]); %[i coef(2),stats(1)]
    indvec=[indvec,i];
    coefvec=[coefvec,coef(2)];
    R2=[R2,stats(1)];
end
[indvec;coefvec;R2]

%%
disp('Table 1')
[~,ind]=sort(vH,'descend');
indnames(indlist(ind(1:10)))
indnames(indlist(ind(end:-1:(end-9))))

disp('Table 2')
[~,ind]=sort(dw,'descend');
indnames(indlist(ind(1:10)))
[~,ind]=sort(etaH,'descend');
indnames(indlist(ind(1:10)))

disp('Table 3: Pearson and Spearman correlations')
corr([etaH',eta',up'])
corr([etaH',eta',up'],'Type','Spearman')

%%
disp('Table OA.1')
[~,ind]=sort(abs(evr(:,1)),'descend');
names=indnames(indlist(ind(1:10)));
[evr_display,ind]=sort(evr(ind(1:10),1),'descend');
names(ind)

[~,ind]=sort(abs(evr(:,2)),'descend');
names=indnames(indlist(ind(1:10)));
[evr_display,ind]=sort(evr(ind(1:10),2),'descend');
names(ind)

[~,ind]=sort(abs(evr(:,3)),'descend');
names=indnames(indlist(ind(1:10)));
[evr_display,ind]=sort(evr(ind(1:10),3),'descend');
names(ind)


%% figure 5
Q=@(t) -b'*sig*(I-sig)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-sig)*t);
tvec=0:.01:5;
qvecSub=zeros(n,length(tvec));
for i=1:length(tvec)
     qvecSub(:,i)=(Q(tvec(i)));
end
%qvecSub=-qvecSub./repmat(qvecSub(:,1),1,length(tvec)); 

qvecSub=qvecSub*100; % this was multiplied by 100 in LT paer
halflives=zeros(n,1);
for i=1:n
    halflives(i)=find(qvecSub(i,:)>-50,1);
end
% one unit of time is a year; multiply the index by 12/100 to convert to months
halflives=halflives*0.12;
%%
try close(fig1); catch; end
ps=100;
fig1=figure(1);
ax=axes; 
hold on; 
xlim([0 max(tvec)]);
ylim([-5 0.1]);

%   {'Construction'                                                   }
%    {'Financial services, except insurance and pension funding'       }
%    {'Computer programming, consultancy and related services'         }
%    {'Employment services'                                            }
%    {'Electricity, transmission and distribution'                     }
%    {'Services auxiliary to financial services and insurance services'}
%    {'Computer, electronic and optical products'                      }
%    {'Real estate services on a fee or contract basis'                }
%    {'Services of head offices; management consulting services'       }
%    {'Advertising and market research services'                       }

%legs = {"Warehousing and support services for transportation" ,"Residential Care and Social Work Activities",...
%"Dyestuffs, agro chemicals", "Accounting, bookkeeping and auditing services; tax consulting services",...
%"Advertising and market research services",...
%"Services auxiliary to financial services and insurance services"}

legs = {"Construction" ,"Financial services, except insurance and pension funding",...
"Computer programming, consultancy and related services",...
"Employment services",...
"Electricity, transmission and distribution",...
"Services auxiliary to financial services and insurance services"}

%
kvec=[57 72 70 86 51 74];

%kvec=[9 10 6 7 8 5];

%line([0,0],[-1,0],'linewidth',1,'linestyle','--','color','black','HandleVisibility','off')
bl=0.2*ones(1,3); gr=0.7*ones(1,3);
plot(tvec,qvecSub(kvec(1),:),'linewidth',2,'linestyle','-','color','#337316');
plot(tvec,qvecSub(kvec(2),:),'linewidth',2,'linestyle',':','color','#055BA6');
plot(tvec,qvecSub(kvec(3),:),'linewidth',2,'linestyle','--','color','#023E73');
plot(tvec,qvecSub(kvec(4),:),'linewidth',2,'linestyle','-','color','#223F59');
plot(tvec,qvecSub(kvec(5),:),'linewidth',2,'linestyle',':','color','#F25922');
plot(tvec,qvecSub(kvec(6),:),'linewidth',2,'linestyle','--','color','#D92344');
%plot([-tvec((ps+1):0:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');
plot([-tvec((ps+1):-1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');

h=legend(legs);
set(h,'FontSize',8,'Location','east');
set(ax,'FontSize', 16, 'xtick',[0 1 2 3 4 5]);
set(ax,'FontSize', 16, 'ytick',[-5 -4 -3 -2 -1 0]);
text(max(tvec)-10,-1.15,"time",'FontSize',18)
xl=xlabel("Years since shock",'FontSize',18,'fontweight','normal');
ylabel(["GDP loss", "(relative to time-zero impact)"],'FontSize',18,'fontweight','normal');
set(gcf,'OuterPosition', [800, 800, 700, 550]);
exportgraphics(fig1,'gdploss_by_sector.pdf');

%%
legs = {"Computer, electronic and optical products","Real estate services on a fee or contract basis"...
"Services of head offices; management consulting services"...
"Advertising and market research services"}

kvec=[39 75 79 82];
try close(fig1); catch; end
ps=100;
fig1=figure(1);ax=axes; hold on; 
xlim([0 max(tvec)]); hold on; 
ylim([-4.5 0.1]);  hold on;
%line([0,0],[-120,0],'linewidth',1,'linestyle','--','color','black','HandleVisibility','off')
bl=0.2*ones(1,3); 
gr=0.7*ones(1,3);
plot(tvec,qvecSub(kvec(1),:),'linewidth',2,'linestyle','-','color','#055BA6');
plot(tvec,qvecSub(kvec(2),:),'linewidth',2,'linestyle',':','color','#023E73');
plot(tvec,qvecSub(kvec(3),:),'linewidth',2,'linestyle','--', 'color','#223F59');
plot(tvec,qvecSub(kvec(4),:),'linewidth',2,'linestyle','-.','color','#D92344');
plot([-tvec((ps+1):-1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');


%plot([tvec((ps+1):1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');
h=legend(legs);
set(h,'FontSize',8,'Location','best');
set(ax,'FontSize', 16, 'xtick',[0 1 2 3 4 5]);
set(ax,'FontSize', 16, 'ytick',[-4 -3 -2 -1 0]);

%ytickformat('percentage');
text(max(tvec)-10,-1.15,"time",'FontSize',18)
xl=xlabel("Years since shock",'FontSize',18,'fontweight','normal');
ylabel(["GDP loss", "(relative to time-zero impact)"],'FontSize',18,'fontweight','normal');
set(gcf,'OuterPosition', [800, 800, 700, 550]);
exportgraphics(fig1,'gdploss_by_sector2.pdf');


%%

kvec=[57 72 70 86 51 74 39 75 79 82];

disp(['half lives'])
mean(halflives)
mean(halflives(kvec(1:3)))
mean(halflives(kvec(4:6)))

disp(['impact after one year'])
mean(qvecSub(kvec(1:3),100))
mean(qvecSub(kvec(4:6),100))

disp(['relative impact after three years'])
mean(qvecSub(kvec(4:6),300))/mean(qvecSub(kvec(1:3),300))

%%  Barrot model extension 
Sig_E = sig;
pfi  = pfi11';

%%  IO matrix
pxij = sig;
pmXi = sum(pxij,2)

%%
% number of sectors
N = length(pxij);
% intermediate good matrix (coefficient)
Omega = diag( pmXi.^(-1) ) * pxij;







%%
%%% eta_i = value added share (of revenue) in sector i
%eta = 1 - pmXi./pyi' ;
eta = 1 - pmXi./pyi ;
    
%%% phi_i = share of final demand i in total output of i
phi = pfi ./ pyi;

%%
%%% Delta_ij =  expenditure on j by i as a share of total production of j
Delta = pxij * diag( pyi.^(-1) ) ;
    
%%% psi_i = share of final demand i in to total final demand
psi = pfi ./sum( pfi);

%%  matrix share of expenduture on j by i in total production j
%normalized the matrix
Gamma=zeros(size(pxij));
for i=1:N
Gamma(i,:) = pxij(i,:) ./ pyi(i);
end

%%  
    
%Labor income share in Value Added 
    
wli_over_vhi = (empcomp./va)';       



%%

    
%% steps to compute the Intermediate good price
% σ is the (constant) elasticity 
% of substitution across goods.
% θ is the (constant) elasticity of substitution between the capital/labor 
% bundle and the intermediate input bundle, 
% and ε is the (constant) elasticity of substitution 
% between intermediate inputs
sigma =  0.9; 
theta =  0.5;
epsilon= 0.001;  
%%

% # Check for zero values before division

int_cons_pyp_nonzero = int_cons_pyp;
int_cons_pyp_nonzer(int_cons_pyp == 0) = 1;% Replace zero with

int_cons_CP_nonzero = int_cons_CP;
int_cons_CP_nonzero(int_cons_CP == 0)= 1 ; % Replace zero with
%

inf_V = int_cons_CP_nonzero./int_cons_pyp_nonzero;

%  assume that:  price_{t-1} was one then price_{t} = inflation_{t}
%



% we can use either the price index of producers
% power of prices
pi_hat = inf_V ;

%pi_hat = [2 3 2 2 3 4 6 2 2 2 3]; % we can use either the price index of producers
%power of prices
pi_hat_1MoinsEpsi = pi_hat.^(1-epsilon);
pi_hat_1MoinsThetai  = pi_hat.^(1-theta);
    
    
%Intermediate Bundle Price
%        Pmi_hat = exp( Omega * log(pi_hat)  );
Pmi_hat = ( Omega .* pi_hat_1MoinsEpsi ).^(1/(1-epsilon));
%%

%%% gamma_i = labor income share in value added in sector i
gamma = wli_over_vhi;

% labor and capital 
li_ai = 0.3 + (1-0.3).*rand(n,1);
ki_ai =  ones(n,1);

%Value added bundle quantity 
ai = (li_ai').^gamma .* (ki_ai').^(1-gamma);
% value added share in total output 
%neta = va/pyi;

%wli_over_GDP
lambda = wli_over_vhi.*va./sum(pfi);
%rli_over_GDP
rho = (1-wli_over_vhi).*va./sum(pfi);
vv = va./ai;


%%   1.95-1

Sub_vec = ones(1, 101);
% Step 2: Set the specified ranges to zeros
Sub_vec(18:20) = 1.95;  % Set rows 18 to 24 to zero
Sub_vec(23:50) = 1.95;  % Set rows 30 to 35 to zero

%%
% Matrix I 
I=eye(n); 


% The parameter d captures the ease of adjustment when input expands
d=0.12; % 12

% data  b is the consumption vector Final consumption expenditure) 
%b = [30 120 90 80 160 80 160 160 300 300 310]'; 
%b = b;  pfi;

b_vec = b; %./sum(b);


shock = +1; % negative or positive
% Update the specification of output
% Q is output (industries intermediate consumption  - intermediate inputs mij ) function 
%Q1=@(t)  ( b_vec' * shock*Sig_E*(I-Sig_E)^-1*expm(-1/d*(I-Sig_E)*t)  )  ;
Q1=@(t)  ( b_vec' * shock*Sig_E*(I-Sub_vec'.*Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t)  )  ;

%% Pmi_hat 
Q2=@(t)  real( ( ( pi_hat_1MoinsEpsi' * (shock *Sig_E*(I-Sub_vec'.*Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t) ) )  ).^(1/(1-epsilon)) );
%% Pi_hat
%Q4=@(t)  real( (  shock^(theta-1) * ( eta'.*ai'.^(1-theta) + (1-eta)'.*( ( pi_hat_1MoinsEpsi'  * (Sig_E*(I-Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t)  )  ).^(1/(1-epsilon))).^(1-theta) ) ).^(1/(1-theta)) );
Q4=@(t)  real( (  shock^(theta-1) * (  eta'.*ai.^(1-theta) +  (1-eta)'.*( ( pi_hat_1MoinsEpsi'  * (Sig_E*(I-Sub_vec'.* Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t)  )  ).^(1/(1-epsilon))).^(1-theta) ) ).^(1/(1-theta)) )
%% pi_hat alterantive (price of final good)
Q5=@(t)  real( (  eta'*(vv).^((theta-1)/theta) + (1-eta)'*( Pmi_hat *shock*Sig_E*(I-Sub_vec'.* Sig_E)^-1*expm(-1./repmat(backlog'/12,n,1).*(I-Sig_E)*t)  ).^((theta-1)/theta)  ).^(theta/(theta-1)) );
%% alternative Pi_hatx
P = 1;
Pi_hat_Alt= ( P.^(1-sigma)./(psi) ).^(1./(1-sigma));

%Q3=@(t)  (  eta*(ai)^((theta-1)/theta) + (1-eta)*( Pmi_hat *shock*Sig_E*(I-Sig_E)^-1*expm(-1/d*(I-Sig_E)*t)  )^((theta-1)/theta)  )^(theta/(theta-1)) ;
%%

% time dimension 20 years
%tvec=0:1:20;  
tvec=0:.01:5;

% vector of sectorl output (3 sectors)
qvecSub1=zeros(n,length(tvec));


for i=1:length(tvec)
     qvecSub1(:,i)=(Q1(tvec(i)));
end
%%

% time dimension 20 years
%tvec=0:1:20;  

% vector of sectorl output (3 sectors)
qvecSub2=zeros(n,length(tvec));


for i=1:length(tvec)
     qvecSub2(:,i)=(Q2(tvec(i)));
end

%%
% vector of sectorl output (3 sectors)
qvecSub4=zeros(n,length(tvec));


for i=1:length(tvec)
     qvecSub4(:,i)=(Q4(tvec(i)));
end

%%
% vector of sectorl output (3 sectors)
qvecSub5=zeros(n,length(tvec));


for i=1:length(tvec)
     qvecSub5(:,i)=(Q5(tvec(i)));
end

%%

%qvecSub1=qvecSub1./repmat(qvecSub1(:,1),1,length(tvec)); 
qvecSub1=qvecSub1*100; % this was multiplied by 100 in LT paer
%%

try close(fig1); catch; end
ps=100;
fig11=figure(1);ax=axes; hold on; xlim([0 max(tvec)]);ylim([-0.1 5]);

legs = {"Construction" ,"Financial services, except insurance and pension funding",...
"Computer programming, consultancy and related services",...
"Employment services",...
"Electricity, transmission and distribution",...
"Services auxiliary to financial services and insurance services"}

%  
kvec=[57 72 70 86 51 74];


%line([0,0],[-1,0],'linewidth',1,'linestyle','--','color','black','HandleVisibility','off')
bl=0.2*ones(1,3); gr=0.7*ones(1,3);
plot(tvec,qvecSub1(kvec(1),:),'linewidth',2,'linestyle','-','color','#337316');
plot(tvec,qvecSub1(kvec(2),:),'linewidth',2,'linestyle',':','color','#055BA6');
plot(tvec,qvecSub1(kvec(3),:),'linewidth',2,'linestyle','--','color','#023E73');
plot(tvec,qvecSub1(kvec(4),:),'linewidth',2,'linestyle','-','color','#223F59');
plot(tvec,qvecSub1(kvec(5),:),'linewidth',2,'linestyle',':','color','#F25922');
plot(tvec,qvecSub1(kvec(6),:),'linewidth',2,'linestyle','--','color','#D92344');
plot([-tvec((ps+1):-1:2),tvec],zeros(1,ps+length(tvec)),'linewidth',1,'color','black');

h=legend(legs);
set(h,'FontSize',8,'Location','east');
set(ax,'FontSize', 16, 'xtick',[0 1 2 3 4 5], 'ytick',[0 1  2  3 4 5]);
%ytickformat('percentage');
text(max(tvec)-10,-1.15,"time",'FontSize',18)
xl=xlabel("Years since shock",'FontSize',18,'fontweight','normal');
ylabel(["GDP Gain", "(relative to time-zero impact)"],'FontSize',18,'fontweight','normal');
set(gcf,'OuterPosition', [800, 800, 700, 550]);
exportgraphics(fig11,'gdpgain_by_sector2SubsisdiesTarget.pdf');


%%


fig1sd = figure('color','w') 
set(fig1sd, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % OR figure('Position', [100, 100, 1024, 1200]);
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 
set(gca,'TickLabelInterpreter','latex')
set(gca,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');



subplot(2,3,1)
plot(tvec,qvecSub1(kvec(1),:),'linewidth',2,'linestyle','-','color','#031CA6'); hold on ;
plot(tvec,qvec1(kvec(1),:),'linewidth',2,'linestyle','-','color','#D92344'); hold on;
title('(a) Construction ','FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

subplot(2,3,2)
plot(tvec,qvecSub1(kvec(3),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(3),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;
title(["(b) Computer programming, ","consultancy and related services "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')





subplot(2,3,3)
plot(tvec,qvecSub1(kvec(2),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(2),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;


title(["(c) Financial services, ","except insurance and pension funding "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')




subplot(2,3,4)
plot(tvec,qvecSub1(kvec(4),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(4),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title('(d) Employment services ','FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

subplot(2,3,5)
plot(tvec,qvecSub1(kvec(5),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(5),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title(["(e) Electricity, transmission"," and distribution "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')




subplot(2,3,6)
plot(tvec,qvecSub1(kvec(6),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(6),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title(["(g) Services auxiliary to financial","services and insurance services"],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

exportgraphics(fig1sd,'gdpsectubsisdiesTarget.pdf');


%%


%%
%try close(fig11); catch; end
%ps=100;
fig22sd = figure('color','w') 
set(fig22sd, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % OR figure('Position', [100, 100, 1024, 1200]);
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 
set(gca,'TickLabelInterpreter','latex')
set(gca,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');



subplot(2,3,1)
plot(tvec,qvecSub2(kvec(1),:),'linewidth',2,'linestyle','-','color','#031CA6'); hold on ;
plot(tvec,qvec2(kvec(1),:),'linewidth',2,'linestyle','-','color','#D92344'); hold on;
title('(a) Construction ','FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Intermediate good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

subplot(2,3,2)
plot(tvec,qvecSub2(kvec(3),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec2(kvec(3),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;
title(["(b) Computer programming, ","consultancy and related services "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Intermediate good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')





subplot(2,3,3)
plot(tvec,qvecSub2(kvec(2),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec2(kvec(2),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;


title(["(c) Financial services, ","except insurance and pension funding "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Intermediate good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')




subplot(2,3,4)
plot(tvec,qvecSub2(kvec(4),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec2(kvec(4),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title('(d) Employment services ','FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Intermediate good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

subplot(2,3,5)
plot(tvec,qvecSub2(kvec(5),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec2(kvec(5),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title(["(e) Electricity, transmission"," and distribution "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Intermediate good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')




subplot(2,3,6)
plot(tvec,qvecSub2(kvec(6),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec2(kvec(6),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title(["(g) Services auxiliary to financial","services and insurance services"],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Intermediate good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')


exportgraphics(fig22sd,'priceintermediatgoodsectubsisdiesTarget.pdf');

%%

hold off;
%%
fig33sd = figure('color','w') 
set(fig33sd, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % OR figure('Position', [100, 100, 1024, 1200]);
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 
set(gca,'TickLabelInterpreter','latex')
set(gca,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');





subplot(2,3,1)
plot(tvec,qvecSub5(kvec(1),:),'linewidth',2,'linestyle','-','color','#031CA6'); hold on ;
plot(tvec,qvec5(kvec(1),:),'linewidth',2,'linestyle','-','color','#D92344'); hold on;
title('(a) Construction ','FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Final good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

subplot(2,3,2)
plot(tvec,qvecSub5(kvec(3),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec5(kvec(3),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;
title(["(b) Computer programming, ","consultancy and related services "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Final good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')





subplot(2,3,3)
plot(tvec,qvecSub5(kvec(2),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec5(kvec(2),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;


title(["(c) Financial services, ","except insurance and pension funding "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Final good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')




subplot(2,3,4)
plot(tvec,qvecSub5(kvec(4),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec5(kvec(4),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title('(d) Employment services ','FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Final good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

subplot(2,3,5)
plot(tvec,qvecSub5(kvec(5),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec5(kvec(5),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title(["(e) Electricity, transmission"," and distribution "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Final good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')




subplot(2,3,6)
plot(tvec,qvecSub5(kvec(6),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec5(kvec(6),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title(["(g) Services auxiliary to financial","services and insurance services"],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Price of Final good", "(relative to time-zero impact)"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')



exportgraphics(fig33sd, 'figalsectubsisdiesTarget.pdf'); % Save as PNG format

hold off
%%  targeted sectors vs target vs no 

kvec=[57 72 70 40 41 42];


       



fig00sd = figure('color','w') 
set(fig00sd, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % OR figure('Position', [100, 100, 1024, 1200]);
set(0,'defaultfigurecolor',[1 1 1]) 
set(0,'defaultaxesfontname','cambria math') 
set(gca,'TickLabelInterpreter','latex')
set(gca,'defaulttextinterpreter','latex')
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');



subplot(2,3,1)
plot(tvec,qvecSub1(kvec(1),:),'linewidth',2,'linestyle','-','color','#031CA6'); hold on ;
plot(tvec,qvec1(kvec(1),:),'linewidth',2,'linestyle','-','color','#D92344'); hold on;
title('(a) Construction ','FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Non-Targeted sectors","output gain"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

subplot(2,3,2)
plot(tvec,qvecSub1(kvec(3),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(3),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;
title(["(b) Computer programming, ","consultancy and related services "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')





subplot(2,3,3)
plot(tvec,qvecSub1(kvec(2),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(2),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;


title(["(c) Financial services, ","except insurance and pension funding "],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
hold on;
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')


%Electrical equipment              
%Machinery and equipment n.e.c.            
%Motor vehicles, trailers and semi-trailers    


subplot(2,3,4)
plot(tvec,qvecSub1(kvec(4),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(4),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title('(d) Electrical equipment','FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["Targeted sectors","output gain"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

subplot(2,3,5)
plot(tvec,qvecSub1(kvec(5),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(5),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title(["(e) Machinery and equipment n.e.c."],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')




subplot(2,3,6)
plot(tvec,qvecSub1(kvec(6),:),'linewidth',2,'linestyle','-','color','#031CA6');  hold on;
plot(tvec,qvec1(kvec(6),:),'linewidth',2,'linestyle','-','color','#D92344');  hold on;

title(["(g) Motor vehicles, ","trailers and semi-trailers"],'FontSize',12,'interpreter','latex')
xl=xlabel("Years since shock",'FontSize',12,'fontweight','normal');
ylabel(["output gain"],'FontSize',12,'fontweight','normal');
legend({'Subsidy program','Baseline'},'FontSize',14,'Location','best','interpreter','latex')

exportgraphics(fig00sd,'gdpsectubsisdiesTarget1.pdf');




%%  table 


indnamesxx=readtable('IO_product_list.xlsx','Sheet', 'Sheet1','ReadVariableNames',false);
indxx=~isnan(table2array(indnamesxx(:,1)));
indnamesxx=table2cell(indnamesxx(indxx,2));

%names_xx = setdiff(indnamesxx, misind, 'stable');



names_xx = indnamesxx;
names_xx(misind) = [];

column2 = qvec1(:,1)

column3 = qvecSub1(:,1)
columndiff = qvecSub1(:,1) - qvec1(:,1)
% Create z1 as a table
z1 = [column2, column3, columndiff];

z1_table = table(z1(:,1), z1(:,2), z1(:,3), 'VariableNames', {'baseline', 'counterfca Scenario', 'Difference'});

% Create z2 as a table or cell array (for names_xx)
z2 = table(names_xx);

% Save to Excel (create a new Excel file or overwrite if exists)
filename = 'output.xlsx';

% Write z1 to the first sheet
writetable(z1_table, filename, 'Sheet', 1);

% Write z2 to the second sheet
writetable(z2, filename, 'Sheet', 2);
