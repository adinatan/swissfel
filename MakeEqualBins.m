function [edges,centres] = MakeEqualBins(x,n,DoPlot)
sortx=sort(x);
edges=zeros((round(length(x)/n))+1,1);
edges(1)=min(x);
edges(end)=max(x);
for ii=2:(round(length(x)/n))
    edges(ii)=sortx((ii-1)*n);
end
if DoPlot
    figure(500)
    hold on
    [bandwidth,density,xmesh,cdf]=kde(x);
    plot(xmesh,density.*n./max(density),'m')
    xlabel('t / ps')
    title({'Histogram'; [num2str(length(edges-1)) ' Bins, with ' num2str(n) ' events in each']})
    for ii=1:numel(edges)
        line([edges(ii) edges(ii)],[0 1.1*n])
    end
    plot(xmesh,density.*n./max(density),'m-','linewidth',1.5)
    axis tight
    legend('kde histogram','bins')
end
centres=edges(1:end-1)+diff(edges)/2;
edges(1)=min(x)-0.1.*(edges(2)-edges(1));% larger bin 1 to keep first point
end

