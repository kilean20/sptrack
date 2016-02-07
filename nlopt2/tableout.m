a=load('table.out');
a1=a(a(:,21)==5000,:);
a2=a(a(:,23)==5000,:);
a3=a(a(:,25)==5000,:);
a4=a(a(:,27)>4500,:);

%%
figure(1)
hold on;
x=[1:26];
plot(x,a(:,21)/5000,'r');
plot(x,a(:,22)/1e-6,'b');
plot(x,a(:,23)/5000,'r');
plot(x,a(:,24)/1e-6,'b');
plot(x,a(:,25)/5000,'r');
plot(x,a(:,26)/1e-6,'b');
plot(x,a(:,27)/5000,'r');
plot(x,a(:,28)/1e-6,'b');
xlabel('cases')
ylabel('emittance growth rate/particle survival rate')
hold off
print( 'emits.eps', '-depsc2')

%%
figure(2)
x=a3(:,26)*1e6;
y=a3;
%x=a4(:,28)*1e6;
%y=a4;
hold on
plot(x,y(:,1) ,'or');
plot(x,y(:,2) ,'og');
plot(x,y(:,3) ,'ob');
plot(x,y(:,4) ,'oc');
plot(x,y(:,5) ,'om');
plot(x,y(:,6) ,'oy');
plot(x,y(:,7) ,'xk');
plot(x,y(:,8) ,'ok');
plot(x,y(:,9) ,'.k');
plot(x,y(:,10),'.k');
h=legend('1','2','3','4','5','6','7','8','9','10');
set(h,'Location','North','Orientation','Horizontal');
xlabel('final emittance (1000 turn) (\mu m)')
ylabel('harmonics')
print( 'harms.eps', '-depsc2')
