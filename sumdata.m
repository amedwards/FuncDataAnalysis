% load \\hscs-share2\centralroot\Research\CAMA\Amanda\FuncDataAnalysis\Dougtable6wk Dougtable6wk
% data=Dougtable6wk(:);
load \\hscs-share2\centralroot\Research\CAMA\Amanda\FuncDataAnalysis\Dougtable12wk Dougtable12wk
nv=length(Dougtable12wk);
columns={'Inst','ID','EGA','BW','Week','Variable','Samples','Week Average'}';
nc=length(columns);
nv=length(Dougtable12wk);
X=zeros(0,nc);
vname=cell(nv,1);
for i=1:length(Dougtable12wk)
    vname{i}=Dougtable12wk(i).variable;
    x=Dougtable12wk(i).data;
    good=~isnan(x(:,5));
    x=x(good,:);
    n=size(x,1);    
    x=[x(:,[1 2 3 4 7]) i*ones(n,1) x(:,[6 5])];
    X=[X;x];
    clear x good
end
g=X(:,1);
id=10*X(:,2)+g;
ega=X(:,3);
bwt=X(:,4);
w=X(:,5);
v=X(:,6);
ns=X(:,7);
y=X(:,8);

ndem=2;
ng=max(g);
nw=max(w);

emax=28;
smin=1000;
pop=ega<=emax&ns>=smin;

sumdata=[];

for k=1:nv

vsub=pop&v==k;
n=zeros(nw,ng);
u=zeros(nw,ng);
s=zeros(nw,ng);
for i=1:nw
    wsub=vsub&w==i;
    for j=1:ng
        gsub=wsub&g==j;
        yy=y(gsub);
        n(i,j)=length(yy);
        u(i,j)=mean(yy);
        s(i,j)=std(yy);
    end
end

sumdata(k).n=n;
sumdata(k).u=u;
sumdata(k).s=s;

end

clear i j k an n u s x good vsub wsub gsub yy

%vplot=[1 3 4 6 7 9 10 12 13 14 15 16];
vplot=[1 3 4 6];
np=length(vplot);

for f=1:np

k=vplot(f);
u=sumdata(k).u;
d=1.96*sumdata(k).s./sqrt(sumdata(k).n);

figure(f)
clf

col=get(gca,'ColorOrder');
for i=1:size(u,2)
    uu=u(:,i);
    dd=d(:,i);
    cc=col(i,:);
    ww=(1:length(uu))';
    h=fill([ww;flipud(ww)],[(uu-dd);flipud(uu+dd)],cc);
    set(h,'facealpha',.25);    
    hold on        
    plot(ww,uu,'Color',cc)
    plot(ww,uu-dd,'Color',cc,'LineStyle',':')        
    plot(ww,uu+dd,'Color',cc,'LineStyle',':')    
end
title(vname{k})
%title([data(k).variable,' p=',num2str(gp(k),'%6.5f')])
    
end
    


