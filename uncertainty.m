avgD_p = sum(D_p)/5;
avgh_p = sum(h_p)/5;
avgD = sum(D)/5;
avgh = sum(h)/5;
uncertain_size = 0.01;
uncertain_mass = 0.01;
uncertain_temp = 0.001;
uncertain_k = 0.001;
c = 377;
lambda =(-4)*m/1000*avgh/1000*c*k/(pi*(avgD/1000)^2*(avgTstable2-avgTstable))*(4*avgh+avgD)/(4*avgh_p+4*avgD_p);
lnd = [4,1/pi,c,m,avgh/1000,-k,avgD^2/1000000,(avgTstable2-avgTstable),(4*avgh_p+avgD_p)/1000,(4*avgh_p+2*avgD_p)/1000];
ech = log(lnd);
uncertainty_lambda_rtv = sqrt((sum(ech))^2*5+(sum(ech)-log(avgD/1000))^2+...
(sum(ech)+4*log(4*avgh_p/1000+avgD_p/1000)-4*log(4*avgh_p/1000+2*avgD_p/1000))^2+...
(sum(ech)+2*log(4*avgh_p/1000+avgD_p/1000)-2*log(4*avgh_p/1000+2*avgD_p/1000))^2);