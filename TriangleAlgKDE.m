function[TriPassFull,TriPass,TriMedian,TriMad]=TriangleAlgKDE(values,DoPlot)
%% add the option for an upper triangle also if a long tail is present.
%values=IPM3(IsData);
%values=EBeam(IsData);
%values=values(IsData);
%values=IPM3;
finitevec=isfinite(values);
values=values(finitevec);
values=values(:)';
nval=sum(values~=0);
%histbins=linspace(min(values(values~=0)),max(values(values~=0)),nval/250);
%[n,bin]=histc(values(values~=0),histbins);
[bandwidth,density,xmesh,cdf]=kde(values);

%n=smooth(n,4);

%% lower cutoff

x1=xmesh(1);
y1=density(1);%mean(n(1:2));
[maxval,maxloc]=(max(density));
x2=xmesh(maxloc);
y2=maxval;
px=xmesh(1:maxloc);
py=density(1:maxloc)';
d = abs((x2-x1)*(y1-py)-(x1-px)*(y2-y1))/sqrt((x2-x1)^2+(y2-y1)^2);

[cutval1,cutloc1]=max(d);
cutval1=xmesh(cutloc1);
%% upper cutoff

x2b=xmesh(end);
y2b=density(end);
%[maxval,maxloc]=(max(n));
x1b=xmesh(maxloc);
y1b=maxval;
pxb=xmesh(maxloc:end);
pyb=density(maxloc:end)';
d = abs((x2b-x1b)*(y1b-pyb)-(x1b-pxb)*(y2b-y1b))/sqrt((x2b-x1b)^2+(y2b-y1b)^2);

[cutval2,cutloc2]=max(d);
cutval2=xmesh(cutloc2+maxloc);
%%
TriPassFull=logical(isfinite(values) & values~=0 & values>cutval1 & values<cutval2);
TriMedian=median(values(TriPassFull));
TriMad=mad(values(TriPassFull)');
TriPass=logical(TriPassFull & values>(TriMedian-TriMad) & values<(TriMedian+TriMad));
%TriMedian
%TriMad

if DoPlot
    figure
    hold on
    plot(xmesh,density,'k-')
    
    
    % figure(777)
    % hold on
    % bar(histbins,n)
    line([x1 x2],[y1 y2])
    plot(xmesh(cutloc1),density(cutloc1),'rx','MarkerSize',15,'linewidth',2.5)
    %plot(histbins(cutloc1),n(cutloc1),'rx','MarkerSize',15,'linewidth',2.5)
    
    %bar(histbins,smooth(n,1))
    line([x1b x2b],[y1b y2b])
    plot(cutval2,density(cutloc2+maxloc),'rx','MarkerSize',15,'linewidth',2.5)
end

% sometimes this kde fucks up when the values are integers
%[bandwidth,density2,xmesh2,cdf]=kde(values(TriPassFull).*(1+rand(1,numel(values(TriPassFull)))/10));
[bandwidth,density2,xmesh2,cdf]=kde(values(TriPassFull));
if DoPlot
    plot(xmesh2,density2.*max(density)./max(density2),'r')
end
[bandwidth,density3,xmesh3,cdf]=kde(values(TriPass));

if DoPlot
    plot(xmesh3,density3.*max(density)./max(density3),'m')

%line([TriMedian TriMedian],[0 max(n)],'linewidth',2.5)
% mad
%line([TriMedian-TriMad TriMedian+TriMad],[0.5*max(n) 0.5*max(n)],'linewidth',2.5)
title({'Histogram';['Total= ' num2str(nval) ' Triangle= ' num2str(sum(TriPassFull)) ' TriangleMedian/Mad= ' num2str(sum(TriPass)) ];[' low= ' num2str(cutval1) ' high= ' num2str(cutval2) ' median= ' num2str(TriMedian) ' mad= ' num2str(TriMad)]})
end
if any(finitevec==0)
    TriPassFulltmp=zeros(size(finitevec));
    TriPassFulltmp(finitevec)=TriPassFull;
    TriPassFull=TriPassFulltmp;
    TriPasstmp=zeros(size(finitevec));
    TriPasstmp(finitevec)=TriPass;
    TriPass=TriPasstmp;
end