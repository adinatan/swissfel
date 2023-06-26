clear all
%close all
h5path='/das/work/p21/p21010/swissfel_reduction/Output/Merged/';
fn=dir([h5path '*_nn20.h5']);
%fn=dir([h5path '*merged.h5']);
tic
% water runs are 91 tp 106
% find which files in fn to load:
runn=0;

    %runnum=[153:162 164:180 182:185];%runn(m);
    %runnum=[187:227 229:243];%runn(m);
    %runnum=[153:177];%runn(m);
    % Comp
   % runnum=[153:162 164:180 182:185]
    % Comp+Dechirp
    %runnum=[187:227 229:265] % something like heat. also in 228
    %runnum=[187:227 229:243 247:265]
    %runnum=[247]% 2
    %runnum=[387 389:396];
    %runnum=[420:423]
    %runnum=[471]
% round jet full power
    runnum=[593:604 606:633 635:639 641:692]
    textprogressbar('load')
    % supercomp
 
for n=1:numel(runnum)
try
    
        runnumstr=num2str(runnum(n),'%0.4i');
        matchC = reshape(strfind({fn.name}, runnumstr), size(fn));
        match(n)  = find(~cellfun('isempty', matchC));
        rn(n)=runnum(n);
catch
rn(n)=NaN;
end
end

    % load files one by one:
    dt=5; %fs


    for n=1:numel(match)
        
        if ~isnan(rn(n))
        textprogressbar(n/numel(match)*100)
        H=Loadh5_Swissfel([fn(match(n)).folder '/' fn(match(n)).name]);

        
        tau=H.timedelays-H.timedelays(1);
        good_time_bins=~isnan(tau);
        %t{n}=tau;
        %S0= H.S0;
        %S2= H.S2;
        %S0a(:,:,n)=H.S0;
        %S2a(:,:,n)=H.S2;

        [S0t,tt]=resample(H.S0(:,good_time_bins)',tau(good_time_bins),1/dt,1,1);
        [S2t,tt]=resample(H.S2(:,good_time_bins)',tau(good_time_bins),1/dt,1,1);
        %S0t(:,:,n)=interp1(tau(1:end-1),H.S0(:,1:1:end-1)',tu,'spline')';
        %S2t(:,:,n)=interp1(tau(1:end-1),H.S2(:,1:1:end-1)',tu,'spline')';
        S0{n}=S0t';
        S2{n}=S2t';
        t{n}=tt;
        sps(n)=double(H.shots_per_step);
    end
    end
%s0(:,size(tmap,2)+1:end)=[];
%s2(:,size(tmap,2)+1:end)=[];

    textprogressbar('done!')
    %%
    if 1
  for n=1:numel(match)
        
        if ~isnan(rn(n))
            tmap(n,1:numel(t{n}))=ones(1,numel(t{n})).*(sps(n));
        end
  end

        s0=zeros(numel(H.q),size(tmap,2));
        s2=zeros(numel(H.q),size(tmap,2));

        for n=1:numel(t)
            s0(:,1:numel(t{n}))=S0{n}(:,1:numel(t{n})).*tmap(n,1:numel(t{n}))+s0(:,1:numel(t{n}));
            s2(:,1:numel(t{n}))=S2{n}(:,1:numel(t{n})).*tmap(n,1:numel(t{n}))+s2(:,1:numel(t{n}));
        end

        s0=s0./sum(tmap);
        s2=s2./sum(tmap);
    end

    %s0=sum(bsxfun(@times,S0t,permute(sps,[3 1 2])),3);
    %s2=sum(bsxfun(@times,S2t,permute(sps,[3 1 2])),3);

    %%
    q=H.q;
    % s0=mean(S0t,3);
    % s2=mean(S2t,3);
   delay=[0:dt:dt*(size(tmap,2)-1)];
    figure('Position',[-2000,500,1000,1000])
    subplot(2,1,1)
    imagesc([0:dt:dt*(size(tmap,2)-1)],H.q,s0)
    %imagesc(tu,q,s0)
    %uimagesc(t{1}(1:end-1),q,S0a(:,1:end-1,4));

    set(gca,'YDir','normal','FontSize',36);
    caxis([-3*mad(s0(:)) 3*mad(s0(:))])
    colormap(flipud(brewermap([],"RdBu")))

    ylabel(['Q (' char(197) '^{-1})'] )
    xlabel('Delay (fs)')
    %xlim([0 600])
    ylim([0.3 4.7])
    title(['\DeltaS0 runs ' num2str(runnum(1)) '-' num2str(runnum(end)) ])
    colorbar

    subplot(2,1,2)
    imagesc([0:dt:dt*(size(tmap,2)-1)],H.q,s2)
    %imagesc(tu,q,s2)
    %uimagesc(t{1}(1:end-1),q,S2a(:,1:end-1,4));
    set(gca,'YDir','normal','FontSize',36);
    caxis([-3*mad(s2(:)) 3*mad(s2(:))])
    colormap(flipud(brewermap([],"RdBu")))
    colorbar
    ylabel(['Q (' char(197) '^{-1})'] )
    xlabel('Delay (fs)')
    ylim([0.3 4.7])
    %xlim([0 600])
    title(['\DeltaS2 runs ' num2str(runnum(1)) '-' num2str(runnum(end)) ])
subtitle(['\sigma=' num2str(estimate_noise(s2(q>0.5 & q<4,delay<100)))])
    %end
%toc

%%
load castner.mat
clear ref_data refi

dc= mean(diff(castner(:,1)))*1000;
ref_data(:,1)=[fliplr(-(dc:dc:dc*500))' ; castner(:,1)*1000];
ref_data(:,2)=[0.*fliplr(-(dc:dc:dc*500))' ; castner(:,2)];
refi=interp1(ref_data(:,1),ref_data(:,2),delay-200,'makima');
plot(delay,refi)

close all
figure('Position',[-2000,500,1000,1000])
%smooth_data=smooth(mean(abs(s2(q>1.6 & q<3.3,:) )),'sgolay',4);
%d=(smooth_data-trimmean(smooth_data(delay<200),66))./...
%    (max(smooth_data(delay<600))-trimmean(smooth_data(delay<200),66));
gfilt=@(x,fwhm)  exp(-((x.^2)/(2*(fwhm/(2*sqrt(2*log(2)))).^2) ));
 
 
    d=smooth(mean(abs(1e6*s2(q>=1 & q<=2,:))),'sgolay',2)'; 

options = optimset('MaxFunEvals',9000,'MaxIter',9000,'TolX',1e-16,'TolFun',1e-16);

tic 
init_vals=[0 1 32 0];
LB = [ -1 0.01  15 -150];
UB = [ 1 10  60  150];


 
[x,fval,exitflag] = fminsearchbnd(@costfun ,init_vals,LB,UB,options,d,refi,delay);
fval*100 
x
t=-100*dt:dt:100*dt;
crefi=x(1)+x(2).*conv(refi,gfilt(t+x(4),x(3)),'same');

plot(delay-200+x(end),d,'-o', delay-200+x(end),crefi,'LineWidth',6,'MarkerSize',20)
xlabel('Delay (fs)')
ylabel('Avg \Delta S_2(1>q>2)')
set(gca,'YDir','normal','FontSize',36);  
    xlim([-200 350])
text(150,1,['\Delta\tau=' num2str(x(3)) ' fs'],'FontSize',53)
title('H_2O + 800nm X-ray scattering in Solution')
subtitle('best optical pump-x-ray probe resolution so far! (25 June 2023)')
legend('Measurment','Theory','FontSize',53)
