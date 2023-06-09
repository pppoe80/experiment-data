avgTstable = sum(d2(:,3))/5;
avgTstable2 = sum(d2(:,2))/5;
time = [30,60,90,120,150,180,210,240,270,300,330,360,390,420,450,480,510,540,570,600,630,660,690,720,750,780,810,840];
T_2 = [65.625,65.125,64.688,64.188,63.688,63.250,62.812,62.375,61.938,61.500,61.062,60.165,60.125,59.188,59.250,58.812,58.375,57.938,57.500,57.062,56.625,56.250,55.938,55.562,55.188,54.572,54.438,54.412];
ts = timeseries(T_2,time);
grid on;hold on;
k = fittedmodel.b*avgTstable;
tavgTstable = log(avgTstable/fittedmodel.a)/fittedmodel.b;
link = @(x)(k*(x-tavgTstable)+avgTstable);
plot(tavgTstable,link(tavgTstable),'o','Color','r');
plot(time,link(time))
plot(time,fittedmodel(time),'Color','b');