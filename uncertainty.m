avgD_p = sum(D_p)/5;
avgh_p = sum(h_p)/5;
avgD = sum(D)/5;
avgh = sum(h)/5;
uncertain_D_a = sqrt(var(D));
uncertain_Dp_a = sqrt(var(D_p));
uncertain_h_a = sqrt(var(h));
uncertain_hp_a = sqrt(var(h_p));
uncertain_size_b = 0.01;
uncertain_mass = 0.01;
uncertain_temp_b = 0.001;
uncertain_temp1 = sqrt(uncertain_temp_b^2+var(d2(:,2)));
uncertain_temp2 = sqrt(uncertain_temp_b^2+var(d2(:,3)));
uncertain_k = 0.001;
uncertain_D = sqrt(uncertain_D_a^2+uncertain_size_b^2);
uncertain_Dp = sqrt(uncertain_Dp_a^2+uncertain_size_b^2);
uncertain_h = sqrt(uncertain_h_a^2+uncertain_size_b^2);
uncertain_hp = sqrt(uncertain_hp_a^2+uncertain_size_b^2);
c = 377;
lambda =(-4)*m/1000*avgh/1000*c*k/(pi*(avgD/1000)^2*(avgTstable2-avgTstable))*(4*avgh+avgD)/(4*avgh_p+4*avgD_p);
lnd = [4,1/pi,c,m,avgh/1000,-k,avgD^2/1000000,(avgTstable2-avgTstable),(4*avgh_p+avgD_p)/1000,(4*avgh_p+2*avgD_p)/1000];
ech = log(lnd);
uncertainty_lambda_rtv = ...
    sqrt(...
    (uncertain_mass/1000)^2+(2*uncertain_D)^2+uncertain_h^2+uncertain_k^2+ ...
    uncertain_temp1^2*(log(avgTstable2-avgTstable))^2+uncertain_temp2^2*(log(avgTstable2-avgTstable))^2+...
    (uncertain_hp/1000)^2*((4*log(4*avgh_p/1000+avgD_p/1000))-4*log(4*avgh_p/1000+2*avgD_p/1000))^2+...
    (uncertain_Dp/1000)^2*((2*log(4*avgh_p/1000+avgD_p/1000))-log(4*avgh_p/1000+2*avgD_p/1000))^2);
uncertainty_lambda = uncertainty_lambda_rtv*lambda;
