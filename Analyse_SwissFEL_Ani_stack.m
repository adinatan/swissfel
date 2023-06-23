function[]=Analyse_SwissFEL_Ani_stack(ScanList)%(Experiment,ScanList,SaveStr,HDF5Path,SavePath,FigurePath,BinSize)
%% Making reduced swissfel files work in the smalldata psana lcls framework
%clear all
addpath(genpath('D:\Dropbox\Dropbox\2306_swissFEL\Analysis'));
%ThisFilename=
%([num2str(HDF5Path) num2str(Experiment,'%s') '_Run' num2str(ScanNumber,'%0.4i') '.h5'])
% ScanNumber=74
% ThisFilename=['F:\SwissFEL\run_' num2str(ScanNumber,'%0.4i') '_reduced_merged.h5']
% 
% 
% HDF5=Loadh5_Swissfel(ThisFilename)
%
ROIRange=[0.8 3.0];
NeighbourSubtraction=1
ChunkSubtraction=0;
climits=[-0.001 0.001]
limits=[-10 20]
BinSize=200
imgoffset=1000
%BinSize=50e-15;
offset=30

%
%ScanList =   [95:99]
reload=true;

for ss=1:numel(ScanList)
    ScanNumber=ScanList(ss);
    disp(['loading ' num2str(ScanNumber)])
    %try % just try loading files
        if reload || ~exist('HDF5','var')
            ThisFilename=['F:\SwissFEL\run_' num2str(ScanNumber,'%0.4i') '_reduced_merged.h5']
            try
                HDF5=Loadh5_Swissfel(ThisFilename)
            catch
                disp('failed to load')
            end
            
            %ThisFilename=([num2str(HDF5Path) num2str(Experiment,'%s') '_Run' num2str(ScanNumber,'%0.4i') '.h5'])
            %            HDF5=LoadlittleData2(ThisFilename)

            cycles=numel(HDF5.scanvar);
            shots=double(HDF5.shots_per_step);
            vec=double(HDF5.mask_ons(:));
            qbins=HDF5.params.qbins;
            phibins=size(HDF5.all_s2d,2);
            %    Enc=double(HDF5.enc.lasDelay(:));

            IPM4=vec;
            IPM5=vec;
            Ebeam=vec;
            fltpos_ps=vec;
            fltpos=vec;
            ampl=vec;
            fwhm=vec;

            TTcor=HDF5.all_tt(:);

            IsData=IPM4>-100;

            HDF5.DropL=logical(HDF5.mask_ofs(:));
            HDF5.DropX=logical(zeros(size(vec)));

            IsData(isnan(HDF5.DropL) | isnan(HDF5.DropX))=false;

            IsDark= HDF5.DropL & HDF5.DropX;
            IsOff=HDF5.DropL & ~HDF5.DropX;
            IsOn= ~HDF5.DropL & ~HDF5.DropX;

            FromScan=ones(size(vec)).*ScanNumber;
            %AllChiPhi=double(HDF5.all_s2d);
            %AllChiPhi=reshape(HDF5.all_s2d,[qbins,phibins,cycles*shots]);

            %AllChiPhi=reshape(HDF5.all_s2d,[qbins,phibins,cycles*shots]);
            %AllChiPhi=permute(AllChiPhi,[3,2,1]);

            AllChiPhi=reshape(HDF5.all_s2d,[qbins,phibins,cycles*shots]);
        

            % make a scanvar for each value
            ScanVec=repmat(HDF5.scanvar,[1,shots])';
            ScanVec=ScanVec(:);

            nominalTime=ScanVec;

            q=HDF5.q;
            nPixels=ones(qbins,1);%double(HDF5.UserDataCfg.epix10k2M.azav__azav_Cake_Npixel(:));

            nChi=ones(1,qbins);
            nChiPhi=ones(1,qbins,phibins);
            %nChi=double(HDF5.UserDataCfg.epix10k2M.azav__azav_Cake_Npixel);
            %nChiPhi=double(HDF5.UserDataCfg.epix10k2M.azav__azav_Cake_norm);

            reload=false;
            VectorOrientation_script
        else % append
            ThisFilename=['F:\SwissFEL\run_' num2str(ScanNumber,'%0.4i') '_reduced_merged.h5']
            try
                HDF5=Loadh5_Swissfel(ThisFilename)
            catch
                disp('failed to load')
            end
            
            %Enc=[Enc;  double(HDF5.enc.lasDelay(:))];
            %AllChiPhi=cat(3,AllChiPhi,double(HDF5.epix10k2M.azav_azav));
             
            cycles=numel(HDF5.scanvar);
            %cycles=[cycles; numel(HDF5.scanvar)];

            shots=[HDF5.shots_per_step];
            %vec=[vec; double(HDF5.mask_ons(:))];
            vec=[double(HDF5.mask_ons(:))];
            qbins=[HDF5.params.qbins];
            phibins=[size(HDF5.all_s2d,2)];
            %    Enc=double(HDF5.enc.lasDelay(:));

            IPM4=[IPM4; vec];
            IPM5=[IPM5; vec];
            Ebeam=[Ebeam; vec];
            fltpos_ps=[fltpos_ps; vec];
            fltpos=[fltpos; vec];
            ampl=[ampl; vec];
            fwhm=[fwhm; vec];

            TTcor=[TTcor; HDF5.all_tt(:)];


            HDF5.DropL=[logical(HDF5.mask_ofs(:))];
            HDF5.DropX=[logical(zeros(size(vec)))];
            IsData=IPM4>100;
           % IsData(isnan(HDF5.DropL) | isnan(HDF5.DropX))=false;

            IsDark=[IsDark'  (HDF5.DropL & HDF5.DropX)'];
            IsOff=[IsOff' (HDF5.DropL & ~HDF5.DropX)'];
            IsOn= [IsOn' (~HDF5.DropL & ~HDF5.DropX)'];

            FromScan=[FromScan; ones(size(vec)).*ScanNumber];
            disp('1')
            %AllChiPhi=double(HDF5.all_s2d);
            % append better
            %AllChiPhi=[AllChiPhi; reshape(HDF5.all_s2d,[qbins,phibins,cycles*shots])];
            
            %AllChiPhiTmp=reshape(HDF5.all_s2d,[qbins,phibins,cycles*shots]);
            %AllChiPhiTmp=permute(AllChiPhiTmp,[3,2,1]);
            %AllChiPhi=[AllChiPhi; AllChiPhiTmp];
            %clear AllChiPhiTmp
            AllChiPhi=cat(3,AllChiPhi,reshape(HDF5.all_s2d,[qbins,phibins,cycles*shots]));
        


            % make a scanvar for each value
            ScanVecTmp=repmat(HDF5.scanvar,[1,shots])';
            ScanVec=[ScanVec; ScanVecTmp(:)];

            nominalTime=ScanVec;

            q=HDF5.q;
            nPixels=ones(qbins,1);%double(HDF5.UserDataCfg.epix10k2M.azav__azav_Cake_Npixel(:));

            nChi=ones(1,qbins);
            nChiPhi=ones(1,qbins,phibins);
            %nChi=double(HDF5.UserDataCfg.epix10k2M.azav__azav_Cake_Npixel);
            %nChiPhi=double(HDF5.UserDataCfg.epix10k2M.azav__azav_Cake_norm);

            reload=false;
            VectorOrientation_script


        %end
    end
end

%%

disp('done loading first file')
clear HDF5


%% generate delays
Tzero=1421440;
RealTimeEnc=(nominalTime+TTcor-Tzero)/1000;% in ps
%RealTimeEnc=RealTimeEnc.*1e-12;
figure;hist(RealTimeEnc,100)
xlabel('delay (s)')
%% extract AllChi params

cab2(1000:8000)
PhiSize=size(AllChiPhi);

%subtract average dark?
%if sum(IsDark)>100
%    AllChiPhi=bsxfun(@minus,AllChiPhi,nanmean(AllChiPhi(:,:,IsDark),3));
%end
%% Make AllChi ~5s, S(q,t) based on S(q,t,phi)
AllChiFromPhi=zeros(PhiSize([1 3]));
AllChiFromPhiNorm=zeros(PhiSize(1),1);

% to minimise ram, done in loop
%cab2(1000:8000)
for ii=1:PhiSize(2) % weigthed average of the phi bins
    % if the Bin is nan, the count should also be nan
    % i am missing checks weather some bins contain nans when they
    % shouldn't, this would mess up normalisation
    tmpChi=bsxfun(@times,squeeze(AllChiPhi(:,ii,:)),nChiPhi(:,ii));
    tmpChi(isnan(tmpChi))=0;
    %AllChiFromPhi=AllChiFromPhi+bsxfun(@times,squeeze(AllChiPhi(:,ii,:)),nChiPhi(:,ii));
    AllChiFromPhi=AllChiFromPhi+tmpChi;
    
    tmpn=nChiPhi(:,ii);
    tmpn(isnan(tmpn))=0;
    AllChiFromPhiNorm=AllChiFromPhiNorm+nChiPhi(:,ii);
   
    
    
end
%%
AllChiFromPhi=bsxfun(@times,AllChiFromPhi,1./AllChiFromPhiNorm);
AllChi=AllChiFromPhi;
clear tmpChi tmpn AllChiFromPhi


%% Params and ROI
PhiSize=size(AllChiPhi);
qROI=isfinite(nanmean(AllChi,2)); % filter nans and infs
%qROI(1)=0;
AllChi=AllChi(qROI,:);
AllChiPhi=AllChiPhi(qROI,:,:);
q=q(qROI);

nimgq=nChi;
nimgphi=nChiPhi;

nimgphi=nimgphi(:,qROI,:);
nimgq=nimgq(qROI);

figure(3000);plot(q,nanmean(AllChi,2)./max(nanmean(AllChi,2)))

hold on
%plot(refTOT(:,1),refTOT(:,2)./max(refTOT(:,2)),'linewidth',1.5)
%legend('data','Solvent Ref')
figure(3001);plot(q,squeeze(nanmean(AllChiPhi,3)))
%% Filter


IntensitySum=double(nansum(AllChi));
IntensitySum=IntensitySum(:);


IPM4=IntensitySum(:);
IPM5=IntensitySum(:);

IsData=IPM4>100 & IPM5>100 & IntensitySum>1e4;

figure;
subplot(3,1,1)
[N,X]=hist(IPM4,1000);
bar(X,N);title('IPM4')
xlim([0 nanmedian(IPM4(IsData))*2])
ylim([0 max(N(2:end))])

subplot(3,1,2)
[N,X]=hist(IPM5,1000);
bar(X,N);title('IPM5')
xlim([0 nanmedian(IPM5(IsData))*2])
ylim([0 max(N(2:end))])

subplot(3,1,3)
[N,X]=hist(IntensitySum,1000);
bar(X,N);title('IntensitySum')
xlim([0 nanmedian(IntensitySum(IsData))*2])
ylim([0 max(N(2:end))])
clear X N;
%% Simple plot binned by scanvec 
ipm4=IPM5; % choose ipm used for filtering. ipm4 is not the same as IPM4, it was IPM4 in run 16
%close all
issimple=0
if issimple
    figure;
    try
        [TriPassFull,TriPass,TriMedian,TriMad,MIN,MAX]=TriangleAlgKDE_limits(ipm4(IsData)./IntensitySum(IsData),1);
        IsDatatmp=(MIN<(ipm4./IntensitySum) & (ipm4./IntensitySum)<MAX);
    catch
        IsDatatmp=IsData;
    end
    figure(1)
    hold on
    plot(IntensitySum(IsDatatmp),ipm4(IsDatatmp),'r.')
    
    % crude filter
    ThisFit=fit(IntensitySum(IsData),ipm4(IsData),'poly1');
    hold on
    ThisConfint=confint(ThisFit);
    %CutOff=4*mean(ThisConfint(:,2));
    CutOff=4*median(ThisConfint(:,2));
    CutOff=10*median(ThisConfint(:,2));
    refline(ThisFit.p1,ThisFit.p2+CutOff)
    refline(ThisFit.p1,ThisFit.p2)
    refline(ThisFit.p1,ThisFit.p2-CutOff)
    
    y2 = polyval([ThisFit.p1 ThisFit.p2+CutOff],IntensitySum);
    y1 = polyval([ThisFit.p1 ThisFit.p2-CutOff],IntensitySum);
    if sum(y1(:)<ipm4(:) & ipm4(:)<y2(:))>numel(IsData)*0.1
        IsData=IsData(:) & y1(:)<ipm4(:) & ipm4(:)<y2(:);
    end
    plot(IntensitySum(IsData), ipm4(IsData),'r.');xlabel('sum');ylabel('IPM')
    title('intensity filterd by linear IPM / intensitysum corr')
    axis tight
end
%% add another IPM cutoff to check gain modes
% remove upper intensities
IsData=IsData & IPM4<(nanmedian(IPM4)+nanstd(IPM4));
sum(IsData)
%% TT filtering
% try
% disp([num2str(sum(IsData)) '  /  ' num2str(numel(IsData)) ' Isdata fraction'])
% %
% fltpos(fltpos(:)<=0 | fltpos(:)>=1024)=0; %
% fltpos(fltpos(:)<=200 | fltpos(:)>=600)=0; %
% ampl(ampl(:)<=0.01 | ampl(:)>0.5)=0; %
% fwhm(fwhm(:)<0 | fwhm(:)>200)=0; %
% 
% 
% fltpos(fltpos(:)<=0 | fltpos(:)>=1024)=0; %
% fltpos(fltpos(:)<=200 | fltpos(:)>=700)=0; %
% ampl(ampl(:)<=0.01 | ampl(:)>1)=0; %
% fwhm(fwhm(:)<0 | fwhm(:)>200)=0; %
% 
% IsTT=fltpos(:)>0 & ampl(:)>0 & fwhm(:)>0; % reset IsTT
% IsTT=IsTT(:);
% % %[TriPassFull2,TriPass,TriMedian,TriMad]=TriangleAlgKDE(fltpos(IsOn(:) & IsTT(:)),0);
% % %[TriPassFull3,TriPass,TriMedian,TriMad]=TriangleAlgKDE(fwhm(IsOn(:) & IsTT(:)),0);
% % %IsTT(IsOn(:) & IsTT(:))=TriPassFull3 & TriPassFull2;
% % 
% % try
% % [TriPassFull1,TriPass,TriMedian,TriMad]=TriangleAlgKDE(fltpos(IsOn(:) & IsTT(:)),1);
% % title('fltpos')
% % xlim([0 1024])
% % [TriPassFull2,TriPass,TriMedian,TriMad]=TriangleAlgKDE(ampl(IsOn(:) & IsTT(:)),1); % error when many around 0
% % title('ampl')
% % %xlim([0 0.1])
% % [TriPassFull3,TriPass,TriMedian,TriMad]=TriangleAlgKDE(fwhm(IsOn(:) & IsTT(:)),1);
% % title('fwhm')
% % xlim([0 350])
% % %IsTT(IsOn(:) & IsTT(:))=TriPassFull3 & TriPassFull2 & TriPassFull1;
% % IsTT(IsOn(:) & IsTT(:))=TriPassFull3 & TriPassFull2;
% % end
% if filterTT
% IsData(IsOn)=IsData(IsOn) & IsTT(IsOn);
% IsData(IsOn(:) & ~IsTT(:))=false;
% end
% % %[TriPassFull4,TriPass,TriMedian,TriMad]=TriangleAlgKDE(ipm4(IsData),1);
% % %[TriPassFull5,TriPass,TriMedian,TriMad]=TriangleAlgKDE(ipm5(IsData),1);
% % %IsData(IsData)= (TriPassFull4 );
% % 
% % IntensitySum=nansum(AllChi)';
% % 
% % %close all
% % figure;plot(IntensitySum(:), ipm4(:),'k.');xlabel('sum');ylabel('IPM')
% % %[TriPassFull,TriPass,TriMedian,TriMad,MIN,MAX]=TriangleAlgKDE_limits(ipm4(IsData)./IntensitySum(IsData),1);
% % %tmp1=(MIN<(ipm4./IntensitySum) & (ipm4./IntensitySum)<MAX);
% % %ThisFit=fit(IntensitySum(tmp1),ipm4(tmp1),'poly1');
% % ThisFit=fit(IntensitySum(IsData),ipm4(IsData),'poly1');
% % hold on
% % ThisConfint=confint(ThisFit);
% % %CutOff=4*mean(ThisConfint(:,2));
% % %CutOff=0.15%
% % CutOff=1*mean(ThisConfint(:,2));
% % %CutOff=1*median(ThisConfint(:,2));
% % refline(ThisFit.p1,ThisFit.p2+CutOff)
% % refline(ThisFit.p1,ThisFit.p2)
% % refline(ThisFit.p1,ThisFit.p2-CutOff)
% % 
% % y2 = polyval([ThisFit.p1 ThisFit.p2+0.5*CutOff],IntensitySum);
% % y1 = polyval([ThisFit.p1 ThisFit.p2-CutOff],IntensitySum);
% % IsData=IsData(:) & y1(:)<ipm4(:) & ipm4(:)<y2(:);
% % plot(IntensitySum(IsData), ipm4(IsData),'r.');xlabel('sum');ylabel('IPM')
% 
% %SaveFigsPNG([FigurePath 'TT_Stat\'],['TT_Stat_' num2str(ScanNumber) '.mat'])
% %close all
% %else % not issimple
%     fltpos(fltpos(:)<=0 | fltpos(:)>=1024)=0; %
%     fltpos(fltpos(:)<=200 | fltpos(:)>=700)=0; %
%     ampl(ampl(:)<=0.01 | ampl(:)>1)=0; %
%     fwhm(fwhm(:)<0 | fwhm(:)>200)=0; %
%     
%     IsTT=fltpos(:)>0 & ampl(:)>0 & fwhm(:)>0; % reset IsTT
%     IsTT=IsTT(:);
% %end
% 
% disp([num2str(sum(IsOn(:) & IsTT(:))) ' TT data'])
% disp([num2str(sum(IsOn(:) & ~IsTT(:))) ' TT discarded for not having TT'])
% disp([num2str(sum(IsData)) '  /  ' num2str(numel(IsData)) ' Isdata fraction'])
% 
% end
IsTT=false(size(IsData));
[Pass,TriPass,TriMedian,TriMad]=TriangleAlgKDE(TTcor(TTcor<800 & TTcor>-300),1)
IsTT(TTcor<800 & TTcor>-300)=Pass;
%tilefigs
%SaveFigsPNG([FigurePath '\Filters\'],['Scan_' num2str(ScanNumber) '_' num2str(SaveStr)],12,8)

%%
%IsData=IsOn|IsOff;
disp(' ')
IsData=IPM4>0;
IsData=IsData(:);
IsDark=IsDark(:);
IsOff=IsOff(:);% | RealTimeEnc<tzero;%in ps
IsOn=IsOn(:);
disp([num2str(sum(IsOn & IsData)) ' IsOn & IsData'])
disp([num2str(sum(IsOff & IsData)) ' IsOff & IsData'])

%ROIRange=[0.8 3.0];
%ROIRange=[0.7 4];
%ROIRange=[2 4];
%ROIRange=[4.2 4.8];
%ROIRange=[4 5];
SaveAngles=q;
ROI=ROIRange(1)<SaveAngles & SaveAngles<ROIRange(2);
% ROI=logical(ones(size(SaveAngles)));
% ROI([1:100 end-249:end])=0;
AllChiScaled=bsxfun(@times,AllChi,1./nanmean(AllChi(ROI,:)));
if ChunkSubtraction
    IsOff=IsOn & IsData & RealTimeEnc<tzero;%in ps
    disp([num2str(sum(IsOff & IsData)) ' IsOff & IsData'])
end
%% qcor, make all phi bins agree for <t0 anisotropy
% fix low q before scaling
ref=squeeze(nanmean(AllChi(:,IsOff&IsData),2));
%ref(q<0.45)=mean(ref(q>0.45 & q<0.60));
figure;plot(q,ref)
for ii=1:size(AllChiPhi,2)
    disp(['qcor bin ' num2str(ii) ' / ' num2str(size(AllChiPhi,2))])
    qcor=squeeze(nanmean(AllChiPhi(:,ii,:),3))./ref;
    AllChiPhi(:,ii,:)=bsxfun(@times,squeeze(AllChiPhi(:,ii,:)),(1./qcor));
    % mask out crappy pixels
end
%%

%% Scale
%Scaler=nanmean(squeeze(AllChi(ROI,:)));
for ii=1:size(AllChiPhi,2)
    disp(['scaling bin ' num2str(ii) ' / ' num2str(size(AllChiPhi,2))])
    AllChiPhi(:,ii,:)=bsxfun(@times,squeeze(AllChiPhi(:,ii,:)),1./nanmean(squeeze(AllChiPhi(ROI,ii,:))));
end
AllDiff=zeros(size(AllChi));
AllDiff=bsxfun(@minus,AllChiScaled,nanmean(AllChiScaled(:,IsOff(:) & IsData(:)),2)); 

%%
figure;plot(q,mean(AllChi(:,IsData),2))
figure;plot(q,squeeze(nanmean(AllChiPhi(:,:,IsData),3)))
%%
IsData=IsData & IPM4>1e4;
AllDiff=zeros(size(AllChi));
AllDiffPhi=zeros(size(AllChiPhi));
if NeighbourSubtraction
    Index=1:numel(IsData); % all elements
    IsGoodOff=IsOff(:) & IsData(:);

    nOff=20;
    OffIndexes=zeros(numel(IsData),nOff); % for anisotropy
    for ii=Index(~IsDark)%
        tmp=IsGoodOff(1:ii-1);%ensure reverse
        OffBefore=(flipud(tmp));
        OffAfter=(IsGoodOff(ii+1:end));
        % make cumsum of both
        
        OffCumSum=zeros(max(length(OffBefore), length(OffAfter)),1);
        OffCumSum(1:length(OffBefore))=OffCumSum(1:length(OffBefore))+OffBefore(:);
        OffCumSum(1:length(OffAfter))=OffCumSum(1:length(OffAfter))+OffAfter;
        OffCumSum=cumsum(OffCumSum);
        
        loc=find(OffCumSum>=nOff,1,'first');
        
        Both=OffCumSum(loc)>nOff;
        
        if Both
            TheseOff=logical(Index(:)>(ii-loc) & Index(:)<=(ii+loc) & IsGoodOff(:) & Index(:)~=ii);
        elseif ~Both
            TheseOff=logical(Index(:)>=(ii-loc) & Index(:)<=(ii+loc) & IsGoodOff(:) & Index(:)~=ii);
        end
        OffIndexes(ii,:)=find(TheseOff)'; % for anisotropy
        %disp(num2str(ii))
        %disp(num2str(find(TheseOff)))
        AllDiff(:,ii)=AllChiScaled(:,ii)-nanmean(AllChiScaled(:,TheseOff),2);
        AllDiffPhi(:,:,ii)=AllChiPhi(:,:,ii)-nanmean(AllChiPhi(:,:,TheseOff),3);
        if mod(ii,1000)==1
            disp([num2str(ii) ' / ' num2str(numel(IsDark)) ' Making Diffs, nearest neighbour'])
        end
    end
elseif ChunkSubtraction    % fast off subtraction
    IsOff=IsOn & IsData & RealTimeEnc<tzero;%in ps
    


    Index=1:numel(IsData); % all elements

    %xpp05816 encoder period ~6000 events
    chunksOn=4e3;%1600; % 1600 period
    chunksOff=6e3;%3200;
    for ii=1:ceil(numel(IsData)./chunksOn)
        disp(['Offsubtracting in chunks of ' num2str(chunksOn) ' , ' num2str(ii)  ' / ' num2str(ceil(numel(IsData)./chunksOn))])
        idx1=0.5.*(ii-1).*chunksOn;
        idx2=idx1+chunksOn;
        
        TheseOn=logical(IsData(:) & idx1<Index(:) & Index(:)<=idx2); % make diffs of all
        
        idx1=0.5.*(ii-1).*chunksOff;
        idx2=idx1+chunksOff;
        
        TheseOff=logical(IsOff(:) & IsData(:) & idx1<Index(:) & Index(:)<=idx2);
        
        AllDiff(:,TheseOn)=bsxfun(@minus,AllChiScaled(:,TheseOn),nanmean(AllChiScaled(:,TheseOff),2));
        AllDiffPhi(:,:,TheseOn)=bsxfun(@minus,AllChiPhi(:,:,TheseOn),nanmean(AllChiPhi(:,:,TheseOff),3));
        
        %AllDiff(:,ii)=AllChiScaled(:,ii)-nanmean(AllChiScaled(:,TheseOff),2);
        %AllDiffPhi(:,:,ii)=AllChiPhi(:,:,ii)-nanmean(AllChiPhi(:,:,TheseOff),3);
    end
else % simpel subtraction
    AllDiff=bsxfun(@minus,AllChiScaled,nanmean(AllChiScaled(:,IsOff&IsData),2));
        AllDiffPhi=bsxfun(@minus,AllChiPhi,nanmean(AllChiPhi(:,:,IsOff&IsData),3));
    disp('Diffs constructed using a simple mean of off and on')
end
disp('Diffs constructed')
Indexes1D=find(IsData);
ScaleMean=mean(AllChi(ROI,:));

%% Plot  Simple
%close all
figure(3000);plot(q,nanmean(AllChi,2)./max(nanmean(AllChi,2)))

hold on
%plot(refTOT(:,1),refTOT(:,2)./max(refTOT(:,2)),'linewidth',1.5)
legend('data','Solvent Ref')

if numel(ScanVec)>1 & numel(unique(ScanVec))>1
    
    ScanVecU=unique(ScanVec);
    try % logspaced timescan
        try
            HDF5.scan.lxt_ttc;
        catch
            HDF5.scan.lxt;
        end
        if sum(~isnan(ScanVec))>0
            AllDelay=zeros(size(AllChi,1),size(ScanVecU,1));
            
            for ii=1:numel(ScanVecU)
                %ThisBin=ScanVec==ScanVecU(ii) & ipm4>0.05;
                ThisBin=ScanVec==ScanVecU(ii) & ipm4>0.1;
                AllDelay(:,ii)=nanmean(AllDiff(:,ThisBin),2);
            end
            
            
            figure(1000+ScanList(1));p=imagesc(ScanVecU,q,AllDelay)
            xlabel('Scan Motor - imagesc spaced')
            title(['Long Time Scans: '  num2str(ScanList(1)) '-' num2str(ScanList(end))])
            set(gca,'Ydir','normal')
            %set(p,'EdgeColor','none')
            caxis([-1 1].*1e-2)
            colormap jet
            print (gcf,'-r300','-dpng','-painters',[FigurePath, 'lxt_ttc', '_', num2str(ScanList(1)+1000), '.png']);
        end
    catch
        if sum(~isnan(ScanVec))>0
            AllDelay=zeros(size(AllChi,1),size(ScanVecU,1));
            
            for ii=1:numel(ScanVecU)
                ThisBin=ScanVec==ScanVecU(ii) & ipm4>0.1 & IsOn;
                %ThisBin=ScanVec==ScanVecU(ii) & ipm4>1 & IsData;
                %AllDelay(:,ii)=nanmean(AllDiff(:,ThisBin),2);
                AllDelay(:,ii)=nanmedian(AllDiff(:,ThisBin),2);
            end
            
            
            figure(1000+ScanList(1));p=pcolor(ScanVecU,q,AllDelay)
            xlabel('Scan Motor')
            set(p,'EdgeColor','none')
            %caxis([-1 1].*1e-2)
            %caxis([-0.005 0.005])
            caxis([climits])
            %caxis([-0.1 0.1])
            %caxis([-0.05 0.05])
            %caxis([-0.01 0.01])
            colormap jet
            colorbar
            title(['Scanvec Binned Scans: '  num2str(ScanList(1)) '-' num2str(ScanList(end))])
        end
    end
    try
        figure(1800+ScanList(1));plot(q,AllDelay)
        legend(num2str(ScanVecU))
        
        %figure(1900+ScanList(1));plot(ScanVecU,nansum(abs(AllDelay(20:end-20,:))))
        simpleroi=q<0.4 & q<1.2;
        simpleroi=q<2 & q<3;
        simpleroi=1.5<q & q<3.2;% water sension
        simpleroi=1.3<q & q<1.7;% ACN x473
        simpleroi=1.<q & q<3;% ACN x473
        figure(1900+ScanList(1));plot(ScanVecU,nansum(abs(AllDelay(simpleroi,:))))
        title('nansum(abs(AllDelay))')
    end
    try %compare
        figure;plot(q,nanmean(AllDelay(),2),'m-','linewidth',2);
        hold on
        acn=load('/cds/home/t/timbvd/Code/XCS-XDS/references/ACN.mat')
        acn=load('/cds/home/t/timbvd/Code/XCS-XDS/references/H2O.mat')
        acn.ACNDiff=acn.H2ODiff;
        %plot(acn.ACNDiff(:,1).*0.83,(acn.ACNDiff(:,2)./max(acn.ACNDiff(acn.ACNDiff(:,1)>1.5,2)))*max(nanmean(AllTTDelayS0(q>1.5&q<3,centres>0),2)),'k--','linewidth',1.5)
        plot(acn.ACNDiff(:,1),(acn.ACNDiff(:,3)./max(acn.ACNDiff(acn.ACNDiff(:,1)>1.5,3)))*max(nanmean(AllDelay,2)),'k--','linewidth',1.5)
        %ylim([climits]*2)
        ylim(climits*2)
    end
% SaveFigsPNG([FigurePath ''],['Simple_' num2str(ScanNumber) '_' num2str(SaveStr)],12,8)   
else
    disp('this is not a simple scan')
end
disp('done plotting simple')

%% XAS
% close all
% 
% %figure(3000);plot(q,nanmean(AllChi,2)./max(nanmean(AllChi,2)))
% 
% %hold on
% %plot(refTOT(:,1),refTOT(:,2)./max(refTOT(:,2)),'linewidth',1.5)
% %legend('data','Solvent Ref')
% 
% if numel(ScanVec)>1 & numel(unique(ScanVec))>1
%     
%     ScanVecU=unique(ScanVec);
%     try % logspaced timescan
%         try
%             HDF5.scan.lxt_ttc;
%         catch
%             HDF5.scan.lxt;
%         end
%         if sum(~isnan(ScanVec))>0
%             AllXASOn=zeros(size(AllChi,1),size(ScanVecU,1));
%             AllXASOff=zeros(size(AllChi,1),size(ScanVecU,1));
%             
%             for ii=1:numel(ScanVecU)
%                 %ThisBin=ScanVec==ScanVecU(ii) & ipm4>0.05;
%                 ThisBin=ScanVec==ScanVecU(ii) & ipm4>0.1;
%                 AllDelay(:,ii)=nanmean(AllDiff(:,ThisBin),2);
%                 
%             end
%             
%             
%             figure(1000+ScanList(1));p=imagesc(ScanVecU,q,AllDelay)
%             xlabel('Scan Motor - imagesc spaced')
%             title(['Long Time Scans: '  num2str(ScanList(1)) '-' num2str(ScanList(end))])
%             set(gca,'Ydir','normal')
%             %set(p,'EdgeColor','none')
%             caxis([-1 1].*1e-2)
%             colormap jet
%             print (gcf,'-r300','-dpng','-painters',[FigurePath, 'lxt_ttc', '_', num2str(ScanList(1)+1000), '.png']);
%         end
%     catch
%         if sum(~isnan(ScanVec))>0
%             AllDelay=zeros(size(AllChi,1),size(ScanVecU,1));
%             
%             for ii=1:numel(ScanVecU)
%                 ThisBin=ScanVec==ScanVecU(ii) & ipm4>0.1 & IsOn;
%                 %ThisBin=ScanVec==ScanVecU(ii) & ipm4>1 & IsData;
%                 %AllDelay(:,ii)=nanmean(AllDiff(:,ThisBin),2);
%                 AllDelay(:,ii)=nanmedian(AllDiff(:,ThisBin),2);
%             end
%             
%             
%             figure(1000+ScanList(1));p=pcolor(ScanVecU,q,AllDelay)
%             xlabel('Scan Motor')
%             set(p,'EdgeColor','none')
%             caxis([-1 1].*1e-2)
%             caxis([-0.005 0.005])
%             %caxis([-0.1 0.1])
%             %caxis([-0.05 0.05])
%             %caxis([-0.01 0.01])
%             colormap jet
%             colorbar
%             title(['Scanvec Binned Scans: '  num2str(ScanList(1)) '-' num2str(ScanList(end))])
%         end
%     end
%     try
%         figure(1800+ScanList(1));plot(q,AllDelay)
%         legend(num2str(ScanVecU))
%         
%         %figure(1900+ScanList(1));plot(ScanVecU,nansum(abs(AllDelay(20:end-20,:))))
%         simpleroi=q<0.4 & q<1.2;
%         simpleroi=q<2 & q<3;
%         simpleroi=1.5<q & q<3.2;% water sension
%         simpleroi=1.3<q & q<1.7;% ACN x473
%         simpleroi=1.<q & q<3;% ACN x473
%         figure(1900+ScanList(1));plot(ScanVecU,nansum(abs(AllDelay(simpleroi,:))))
%         title('nansum(abs(AllDelay))')
%     end
%     try %compare
%         figure;plot(q,nanmean(AllDelay(),2),'m-','linewidth',2);
%         hold on
%         acn=load('/cds/home/t/timbvd/Code/XCS-XDS/references/ACN.mat')
%         acn=load('/cds/home/t/timbvd/Code/XCS-XDS/references/H2O.mat')
%         acn.ACNDiff=acn.H2ODiff;
%         %plot(acn.ACNDiff(:,1).*0.83,(acn.ACNDiff(:,2)./max(acn.ACNDiff(acn.ACNDiff(:,1)>1.5,2)))*max(nanmean(AllTTDelayS0(q>1.5&q<3,centres>0),2)),'k--','linewidth',1.5)
%         plot(acn.ACNDiff(:,1),(acn.ACNDiff(:,3)./max(acn.ACNDiff(acn.ACNDiff(:,1)>1.5,3)))*max(nanmean(AllDelay,2)),'k--','linewidth',1.5)
%         %ylim([climits]*2)
%         ylim(climits*2)
%     end
% % SaveFigsPNG([FigurePath ''],['Simple_' num2str(ScanNumber) '_' num2str(SaveStr)],12,8)   
% else
%     disp('this is not a simple scan')
% end
% disp('done plotting simple')

%% Plot TT

%close all
Norm=IntensitySum;
%BinSize=0.020e-12;
%BinSize=500;

IsTT=IsData
%cc
%IpmCut=cutlist(cc);
IsE=IsData;
sum(IsE)

if (sum(IsTT)./numel(IsTT))>0.25
    if BinSize<1
        [edges,centres] = MakeLinearBins(RealTimeEnc(IsData & IsOn & IsE & RealTimeEnc>limits(1) & RealTimeEnc<limits(2)),BinSize,0);
    else
        %BinSize=200;
        
        [edges,centres] = MakeEqualBins(RealTimeEnc(IsData & IsOn & IsE & RealTimeEnc>limits(1) & RealTimeEnc<limits(2)),BinSize,1);disp([num2str(sum(IsOff & IsData)) ' IsOff & IsData'])
        if numel(edges)>numel(unique(edges))
            edges=unique(edges);
            centres=edges(1:end-1)+(diff(edges)/2);
        end
    end
    %[edges,centres]=MakeLinearBins([-2 2],0.010,0);
    %[edges,centres]=MakeLinearBins([-1.2 -0.4],0.025,0);
    %[edges,centres]=MakeLinearBins([-50 50],1,0);
    %[edges,centres]=MakeLinearBins([-2 6],0.1,0);
    else
    centres=ScanVecU;
    edges=ScanVecU-100e-15;
    edges(end+1)=edges(end)+200e-15    
    RealTimeEnc=ScanVec;
end
newhsv=flipud(cool(numel(centres)));
AllTTDelay=zeros(size(AllChi,1),size(centres,1));
%close all
%AllTTXAS=zeros(2,size(centres,1));
%AllTTXASOff=zeros(2,size(centres,1));
%DataIOff=zeros(size(centres));
%ErrI=zeros(size(centres));
%ErrIOff=zeros(size(centres));
nImg=zeros(numel(centres),2);
figure(200)
hold on
for ii=1:numel(centres)
    BinOn=edges(ii)<RealTimeEnc & edges(ii+1)>=RealTimeEnc & IsData & IsOn & IsE;
    if sum(IsOff)>10
        BinOff=edges(ii)<RealTimeEnc & edges(ii+1)>=RealTimeEnc & IsData & IsOff & IsE;
        %DataIOff(ii)=mean(sum(Dio(BinOff))./sum(Norm(BinOff)));
        %ErrIOff(ii)=std(((Dio(BinOff))./(Norm(BinOff))));
        %AllTTXASOff(:,ii)=[nansum(DiodeU2(BinOff))./nansum(IntensitySum(BinOff)) nansum(DiodeU4(BinOff))/nansum(IntensitySum(BinOff))];
        %AllTTXASOff(:,ii)=[nansum(DiodeU2(BinOff))./nansum(ipm4(BinOff)) nansum(DiodeU4(BinOff))/nansum(ipm4(BinOff))];
    else
        BinOff=[];
    end
    AllTTDelay(:,ii)=nanmean(AllDiff(:,BinOn),2);
    %AllTTDelay(:,ii)=nanmean(AllDiff(:,BinOn),2)./sum(Norm(BinOn));
    %ErrI(ii)=std(((Dio(BinOn))./(Norm(BinOn))));
    nImg(ii,:)=[sum(BinOn) sum(BinOff)];
    plot(q,AllTTDelay(:,ii)*100,'color',newhsv(ii,:),'linewidth',1.5)
    %AllTTXAS(:,ii)=[nansum(DiodeU2(BinOn))./nansum(IntensitySum(BinOn)) nansum(DiodeU4(BinOn))/nansum(IntensitySum(BinOn))];
    %AllTTXAS(:,ii)=[nansum(DiodeU2(BinOn))./nansum(ipm4(BinOn)) nansum(DiodeU4(BinOn))/nansum(ipm4(BinOn))];
    %axis tight
    %drawnow
    %ylim
end
axis tight
ylim(climits.*500)
legend(num2str(centres))
figure;
subplot(3,1,1:2)
h=pcolor(centres,q,AllTTDelay)
set(h,'EdgeColor','none')
caxis([-0.5 0.5].*1e-2)
xlabel('Delay (ps)')
ylabel('Q (A-1)')
hold on
ylimits=get(gca,'ylim')
%plot([tzero tzero],ylimits,'b-')
title(['Run ' num2str(ScanList(:)')])

subplot(3,1,3)
bar(centres, nImg(:,1),'r')
hold on
bar(centres, nImg(:,2),'k')
ylimits=get(gca,'ylim')
%plot([tzero tzero],ylimits,'b-')
legend('on','off')
title('n events')
 %SaveFigsPNG([FigurePath ''],['SimpleTT_' num2str(ScanNumber) '_' num2str(SaveStr)],12,8)  
% %%
kinroi=[0.5 2.5]
figure(2000+ScanList(1))
plot(centres,nanmean(abs(AllTTDelay(kinroi(1)<q & q<kinroi(2),:)),1)./max(nanmean(abs(AllTTDelay(kinroi(1)<q & q<kinroi(2),:)),1)))
hold on
plot(centres,medfilt1(nanmean(abs(AllTTDelay(kinroi(1)<q & q<kinroi(1),:)),1),1),'linewidth',1.5)
ylimits=get(gca,'ylim')
%plot([tzero tzero],ylimits,'r-')
title(['Run ' num2str(ScanList(:)')])
% add svd kin
%hold on



%SaveFigsPNG([FigurePath ''],['SimpleTT_' num2str(ScanNumber) '_' num2str(SaveStr)],12,8)  
Delays=centres;
%save([SavePath 'Scan_' num2str(ScanNumber) '.mat'],'Delays','q','AllTTDelay','nImg','ScanList')
a=medfilt1(nanmean(abs(AllTTDelay(kinroi(1)<q & q<kinroi(2),:)),1),5);
refline(0,nanmean(a(centres<1)))
xlabel('Delay (ps)')

tilefigs
%% make binned AllDiffPhi
AllPhiDelay=zeros(size(AllChi,1),size(AllChiPhi,2),size(centres,1));
for ii=1:numel(centres)
    BinOn=edges(ii)<RealTimeEnc & edges(ii+1)>=RealTimeEnc & IsData & IsOn & IsE;
    if sum(IsOff)>10
        BinOff=edges(ii)<RealTimeEnc & edges(ii+1)>=RealTimeEnc & IsData & IsOff & IsE;
    else
        BinOff=0;
    end
    AllPhiDelay(:,:,ii)=nanmean(AllDiffPhi(:,:,BinOn),3);
    disp([num2str(ii) ' / ' num2str(numel(centres))])
end


%% phibins for anisotropy
%climits=climits/12
%SaveStr='auto'
%close all
cab2(1000:8000)
% scale and subtract off
% Lorentz seperation
AllTTDelayS0=zeros(size(AllTTDelay)).*nan;
AllTTDelayS2=zeros(size(AllTTDelay)).*nan;
nPhi=size(AllDiffPhi,2);



a=([0:(nPhi-1)].*(360/nPhi))+offset;
a(a>180)=-360+a(a>180);
%r = Q2TT(q,8.8);
r = Q2TT(q,9.5);
nbinR=numel(q);

for tt=1:numel(centres)
    disp(['Seperationg S2/S0 for bin ' num2str(tt) ' / ' num2str(numel(centres))])
    S=zeros(nbinR,2);
    for qq=1:numel(q)
                ThisI=AllPhiDelay(qq,:,tt);
        % r=q
        % a=azimuth
        Thisr=r(qq);%q()
        
        if sum(~isnan(ThisI) & ThisI~=0)>3
            PhiROI=~isnan(ThisI) & ThisI~=0;
            Thisa=a(PhiROI);
            ThisI=ThisI(PhiROI);
            
            X=-cosd(Thisr/2).*cosd(Thisa);
            %X=-cosd(Thisa);
            P2=0.5*((3*X.^2)-1);
            

            [A B Stat]=Theil_Sen_Regress_Stat(P2',ThisI');
            %line([min(P2) max(P2)],[A.*min(P2)+B  A.*max(P2)+B]);
            %[A,B]=weightedfit([P2' ThisI' (1./sqrt(nimgphi(qq,PhiROI)))']);
            %B=nanmean(ThisI);
            %A=Results.slope;
            %B=Results.Intercept;
            
            
            %Stats(qq,:)=Stat;
            S(qq,:)=[A,B];
%             plot(P2,ThisI,'b.')
%             
            % plot
%             figure(1);plot(P2,ThisI,'k*')
%             hold on
%             line([min(P2) max(P2)],[A.*min(P2)+B  A.*max(P2)+B],'color',[1 0 0]);
%             title(num2str(q(qq)))
%             hold off
%             pause(0.1)
            
            %xlabel('P2(cos(\theta_Q))')
            %ylabel('\DeltaS (a.u.)')
        end
    end
    AllTTDelayS0(:,tt)=S(:,2);
    AllTTDelayS2(:,tt)=S(:,1);
end

%close all

disp([num2str(tt) ' / ' num2str(numel(centres))])


    figure(500+ScanNumber+imgoffset);h=pcolor(centres,q,AllTTDelay)
    set(h,'EdgeColor','none')
    %caxis([climits].*1e-2)
    caxis([climits])
    title(['Integrated, scan ' num2str(ScanList(1)) '-' num2str(ScanList(end))])
    xlabel('Delay (ps)')
    ylabel('Q (A^-^1)')
    colormap jet
    colorbar
    
    figure(2*ScanNumber+1000+imgoffset);h=pcolor(centres,q,AllTTDelayS0)
    set(h,'EdgeColor','none')
    %caxis([climits].*1e-2)
    caxis([climits])
    title(['S0, scan '  num2str(ScanList(1)) '-' num2str(ScanList(end)) ' ' num2str(sum(IsData&IsOn)) ' / ' num2str(numel(IsData)) ])
    xlabel('Delay (ps)')
    ylabel('Q (A^-^1)')
    colormap jet
    colorbar
    
    figure(2*ScanNumber+1001+imgoffset);h=pcolor(centres,q,AllTTDelayS2)
    set(h,'EdgeColor','none')
    %caxis([climits].*1e-2./2)
    caxis([climits])
    title(['S2, scan '  num2str(ScanList(1)) '-' num2str(ScanList(end)) ' ' num2str(sum(IsData&IsOn)) ' / ' num2str(numel(IsData)) ])
    xlabel('Delay (ps)')
    ylabel('Q (A^-^1)')
    colormap jet
    colorbar
   %save([SavePath 'Scan_' num2str(ScanNumber) '.mat'],'Delays','q','AllTTDelay','centres','AllTTDelayS0','AllTTDelayS2','nImg','ScanList') 
    %%
%     QROI=(q>1 & q<3);
%     QROI=(q>0.8 & q<2.5);
%     figure(ScanNumber+6000+imgoffset);plot(centres,nanmean(AllTTDelayS2(QROI,:)));title(['Anisotropy kinetics, run ' num2str(ScanNumber)])
%     hold on
%     % calculate COMs to detecto oscillations
%     COMs=zeros(size(centres));
%     for ii=1:numel(COMs)
%         Thisx=q(QROI);
%         Thisy=abs(AllTTDelayS2(QROI,ii));
%         COMs(ii)=mean(Thisx(:).*Thisy(:))/sum(Thisy(:));
%     end
%     plot(centres,COMs,'r-')
%     legend('nanmean','COMs')
%     
%     caxis([-0.001 0.001]*2)
%     figure(9999);plot(centres,sum(abs(AllTTDelayS2(q>1&q<3,:)))./max(sum(abs(AllTTDelayS2(q>1&q<3,:)))))
%     
%     %[U,S,V0]=svds(AllTTDelayS0,1);
%     %[U,S,V2]=svds(AllTTDelayS2,1);
%     [U,S,V0]=svds(AllTTDelayS0(QROI,:),1);
%     [U,S,V2]=svds(AllTTDelayS2(QROI,:),1);
%     hold on 
%     plot(centres,V0(:,1),'linewidth',1.5)
%     plot(centres,V2(:,1),'linewidth',1.5)
% %
% %%
% [U,S,V]=svds(AllTTDelay(q>1 & q< 3,2:end-1),1)
% figure;plot(centres(2:end-1),V(:,1))
% figure;plot(q(q>1 & q< 3),U(:,1)./max(U(:,1)))
% hold on
% acn=load('/cds/home/t/timbvd/Code/XCS-XDS/references/H2O.mat')
% acn.ACNDiff=acn.H2ODiff;
% %plot(acn.ACNDiff(:,1).*0.83,(acn.ACNDiff(:,2)./max(acn.ACNDiff(acn.ACNDiff(:,1)>1.5,2)))*max(nanmean(AllTTDelayS0(q>1.5&q<3,centres>0),2)),'k--','linewidth',1.5)
% plot(acn.ACNDiff(:,1),(acn.ACNDiff(:,3)./max(acn.ACNDiff(acn.ACNDiff(:,1)>1.5,3))),'m--','linewidth',1.5)