
Ang = char(197);
%%
% anisotropy from q,phi
% 1. find phi=0 (32 deg?)
 
qsp=q>1.35 & q<1.55;
as=nanmean(nanmean(davg(:,:,qsp),3),1);

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
%plot(phi_grid,as,phi_grid,modelfun(b,phi_grid)); hold on;

% 
phi_corrected=phi_grid+b(end);
phi_corrected=[phi_corrected phi_corrected+2*pi];
ase=[as as];
%%

b_params=[2];
clear betas ese Betas ESE
ang_grid=linspace(0,2*pi,12);
for nt=1:size(davg,1)
for nq=1:size(davg,3)
    davgc(nt,:,nq)=interp1(phi_corrected,    [davg(nt,:,nq) davg(nt,:,nq)], ang_grid,'spline');
    dsavgc(nt,:,nq)=interp1(phi_corrected,    [dsavg(nt,:,nq) dsavg(nt,:,nq)], ang_grid,'spline');

    w=1./dsavgc(nt,:,nq).^2;
    vec= davgc(nt,:,nq);
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

%%
XrayEnergy_in_keV = 16.01; 
% The length of the wave vector of the incoming beam in A^(-1)
k0 = XrayEnergy_in_keV*1e3 / 27.2113966 / 137.036 / 0.5291772;


%sf=@(n) cos ( asin(xq(:)./(2*k0)) ) .^((n-1)*2);%sqrt(1-(qgrid./(2.*k0)).^2);

sf=@(n) cos ( asin(q(:)./(2*k0)) ) .^((n-1)*2);%sqrt(1-(qgrid./(2.*k0)).^2);


clear  sn
for n=1:size(Betas,3)
    sn(:,:,n) = bsxfun(@rdivide, Betas(:,:,1).* Betas(:,:,n),sf(n));
    sn(:,:,1)= Betas(:,:,1);
end
