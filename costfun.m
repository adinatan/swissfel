function cost=costfun(x,d,refi,delay)

dt=mean(diff(delay));
t=-100*dt:dt:100*dt;

filt=@(t,fwhm)  exp(-((t.^2)/(2*(fwhm/(2*sqrt(2*log(2)))).^2) ));
 
crefi=x(1)+x(2).*conv(refi,filt(t+x(4),x(3)),'same');
   
range_fit=(delay>0 & delay<550)  ;
cost= nansum( abs(  d(range_fit)-crefi(range_fit) ) ) ./sum(range_fit);


end
