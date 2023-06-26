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
    laser_power=[4,10,17,24,35];

    totrunnum=reshape(436:459,4,[]);


for m=1:numel(laser_power)
    textprogressbar('load')
runnum=totrunnum(:,m);

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

       s0=zeros(500,100);
       s2=zeros(500,100);

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

        tmap(n,1:numel(t{n}))=ones(1,numel(t{n})).*(sps(n));

        s0(:,1:numel(t{n}))=S0{n}(:,1:numel(t{n})).*tmap(n,1:numel(t{n}))+s0(:,1:numel(t{n}));
        s2(:,1:numel(t{n}))=S2{n}(:,1:numel(t{n})).*tmap(n,1:numel(t{n}))+s2(:,1:numel(t{n}));
     
         end
    end
s0(:,size(tmap,2)+1:end)=[];
s2(:,size(tmap,2)+1:end)=[];


    textprogressbar('done!')
   
        s0=s0./sum(tmap);
        s2=s2./sum(tmap);
 
 vec(m)=mean(mean(abs(s0(H.q>0,:))));
 end
    %%
    q=H.q;
    % s0=mean(S0t,3);
    % s2=mean(S2t,3);
   
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
delay=[0:dt:dt*(size(tmap,2)-1)];
figure('Position',[-2000,500,1000,1000])
plot(delay,mean(abs(s2(q>1.6 & q<2.4,:))),'-s','LineWidth',7,'MarkerSize',22)
  xlabel('Delay (fs)')
 set(gca,'YDir','normal','FontSize',36);  
 %xlim([0 600])
