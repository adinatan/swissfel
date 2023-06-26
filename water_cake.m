clear all
close all
h5path='/das/work/p21/p21010/swissfel_reduction/Output/Merged/';
fn=dir([h5path '*.h5']);
Ang = char(197);
tic
% water runs are 91 tp 106
% find which files in fn to load:
runn=0;

    %runnum=[153:162 164:171];%runn(m);
    runnum=[153:162 164:178:180 182:185];%runn(m);
    %runnum=[187:225];%runn(m);
    for n=1:numel(runnum)
        runnumstr=num2str(runnum(n),'%0.4i');
        matchC = reshape(strfind({fn.name}, runnumstr), size(fn));
        match(n) = find(~cellfun('isempty', matchC));
    end

    % load files one by one:
    dt=10; %fs


    for n=1:numel(match)
        n/numel(match)*100
        H=Loadh5_Swissfel([fn(match(n)).folder '/' fn(match(n)).name]);

        tau=H.timedelays-H.timedelays(1);
        good_time_bins=~isnan(tau);
        

        %%
        [caket,tt]=resample(permute(H.dS_caked(:,:,good_time_bins),[3 1 2]),tau(good_time_bins),1/dt,1,1);
        %[S2t,tt]=resample(H.S2(:,good_time_bins)',tau(good_time_bins),1/dt,1,1);
        %S0t(:,:,n)=interp1(tau(1:end-1),H.S0(:,1:1:end-1)',tu,'spline')';
        %S2t(:,:,n)=interp1(tau(1:end-1),H.S2(:,1:1:end-1)',tu,'spline')';
        cake{n}=permute(caket,[2 1 3]);
        %S2{n}=S2t';
        t{n}=tt;
        sps(n)=double(H.shots_per_step);
    end

    %%
    if 1
        for n=1:numel(match)
            tmap(n,1:numel(t{n}))=ones(1,numel(t{n})).*(sps(n));
        end

        cakei=zeros(numel(H.q),size(tmap,2),size(caket,3));

        for n=1:numel(t)
            cakei(:,1:numel(t{n}),:)=cake{n}(:,1:numel(t{n}),:).*tmap(n,1:numel(t{n}),:)+cakei(:,1:numel(t{n}),:);
        end

        cakei=cakei./sum(tmap);
        
    end

%%
%%
% anisotropy from q,phi
% 
  
q=H.q;
qsp=q>1.6 & q<2.2;
delay=[0:dt:dt*(size(tmap,2)-1)];%50:60;
ti=delay>100 & delay<500;

as=squeeze(nanmean(nanmean(cakei(qsp,ti,:),2),1));

phi_grid = linspace(0, 2*pi, numel(as)+1);
dphi=mean(diff(phi_grid));
phi_grid = phi_grid(1:end-1) ;

B0 = mean(as);  % Vertical shift
B1 = (max(as) - min(as))/2; % Amplitude
B2 = 2; % Phase (Number of peaks)
B3 = 0; % Pnhase shift (eyeball the Curve)
myFit = NonLinearModel.fit(phi_grid,as, 'y ~ b0 + b1*sin(b2*x1 + b3)', [B0, B1, B2, B3])

modelfun = @(b,x) b(1) + b(2)*sin(x*b(3) +b(4));
b= myFit.Coefficients{:,1};
 plot(phi_grid,as,phi_grid,modelfun(b,phi_grid)); hold on;

% 
phi_corrected=phi_grid+b(end);
phi_corrected=[phi_corrected phi_corrected+2*pi];
ase=[as as];

b_params=[2];
clear betas ese Betas ESE
ang_grid=linspace(0,2*pi,12);
davg=permute(cakei,[2 3 1]);
for nt=1:size(davg,1)
for nq=1:size(davg,3)
    davgc(nt,:,nq)=interp1(phi_corrected,    [davg(nt,:,nq) davg(nt,:,nq)], ang_grid,'makima');
   % dsavgc(nt,:,nq)=interp1(phi_corrected,    [dsavg(nt,:,nq) dsavg(nt,:,nq)], ang_grid,'spline');

    
    vec= davgc(nt,:,nq);
    w=ones(size(vec));%./dsavgc(nt,:,nq).^2;
     if ~all(vec==0) | ~any(isinf(w))
    [betas ]=LDSDw(vec,b_params, w );
     else
betas =NaN(numel(b_params)+1,1);
% ese =[NaN;NaN:NaN];
     end

     Betas(nq,nt,:)=betas;
    % ESE(nq,nt,:)=ese;
     end
end


XrayEnergy_in_keV = 9.0; 
% The length of the wave vector of the incoming beam in A^(-1)
k0 = XrayEnergy_in_keV*1e3 / 27.2113966 / 137.036 / 0.5291772;
sf=@(n) cos ( asin(q(:)./(2*k0)) ) .^((n-1)*2);%sqrt(1-(qgrid./(2.*k0)).^2);


clear  sn
for n=1:size(Betas,3)
    sn(:,:,n) = bsxfun(@rdivide, Betas(:,:,1).* Betas(:,:,n),sf(n));
    sn(:,:,1)= Betas(:,:,1);
end


%%

 
    figure('Position',[-2000,500,1000,1000])
    subplot(2,1,1)
        s0=sn(:,:,1);
    imagesc(delay,q,s0)
    %imagesc(tu,q,s0)
    %uimagesc(t{1}(1:end-1),q,S0a(:,1:end-1,4));

    set(gca,'YDir','normal','FontSize',36);

    caxis([-2*mad(s0(:)) 2*mad(s0(:))])
    colormap(flipud(brewermap([],"RdBu")))

    ylabel(['Q (' char(197) '^{-1})'] )
    xlabel('Delay (fs)')
    %xlim([-200 1200])
    ylim([0.3 4.7])
    title(['\DeltaS0 runs ' num2str(runnum(1)) '-' num2str(runnum(end)) ])
    colorbar

    subplot(2,1,2)
    s2=squeeze(sn(:,:,2));
    imagesc(delay,q,s2)
    %imagesc(tu,q,s2)
    %uimagesc(t{1}(1:end-1),q,S2a(:,1:end-1,4));
    set(gca,'YDir','normal','FontSize',36);
    caxis([-2*mad(s2(:)) 2*mad(s2(:))])
    colormap(flipud(brewermap([],"RdBu")))
    colorbar
    ylabel(['Q (' char(197) '^{-1})'] )
    xlabel('Delay (fs)')
    ylim([0.3 4.7])
    %xlim([-200 1200])
    title(['\DeltaS2 runs ' num2str(runnum(1)) '-' num2str(runnum(end)) ])
%end
%toc

%%
figure
delay=[0:dt:dt*(size(tmap,2)-1)];
figure('Position',[-2000,500,1000,1000])
plot(delay,(mean(abs(s2(q>1.6 & q<2.2,:)))),'-x','LineWidth',7,'MarkerSize',20)
 xlabel('Delay (fs)')
 set(gca,'YDir','normal','FontSize',36);  
xlim([0 600])
title(['Runs ' num2str(runnum(1)) '-' num2str(runnum(end)) ' mean q>1.6 & q<2.2'])









