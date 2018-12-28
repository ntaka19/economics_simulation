clear;
clf;

[x,y] = ndgrid(0:0.01:0.5,0:0.01:3);


%%PARAMETERS
A = 0.15; %A=a*h*l*zeta;
eps = 0.20; %ε
c0 = 0.001;  
c1 = 0.1; %keijo-keihiritu 
c2 = 0.4; 

B = 0.15; %B=a*d*m*zeta
alpha = 6.0; %α = bm 
f0 = -0.59;


%%INTERIOR PARAMETERS
a = 0.625;
hl = 0.3; %h*l
%m = 3;

b = 2.0;
zeta = 0.8;

rho = 0.1;

%%ORIGINAL VALUE
d = 0.1; %hl = md =0.3
mu = 0.60; %μ

%%%%%%%%%%%%%%%%%%%%%%%%g=omega=0
i = B.*y.^eps+b.*m.*x.^mu+f0; 
s = A.*y.^eps-c1.*y-c2.*x-c0;
h = A.*y.^eps-c1.*y-c2.*x-c0-(y+x).*(B.*y.^eps+alpha.*x.^mu+f0);

figure(1);
I = contour(x,y,i,[0,0],'b'); hold on;
%clabel()
S = contour(x,y,s,[0,0],'r'); hold on;
H = contour(x,y,h,[0,0],'g'); hold on;

xlabel('一人あたりの公共生活施設量(ω)');
ylabel('一人あたりの公共生産施設量(g)');


%%%%%%%% DIFFERENTIAL EQUATION %%%%%%%

DT = 0.001;
TT = 10000;
TOUT = 5;

%%INITIAL STATE
g = 0.3;
omega = 0.2;

p1 = 0.3;
p2 = 0.1;


tt(1) = 0.0;
gg(1) = g;
omegalist(1) = omega;
p1list(1) = p1;
p2list(1) = p2;


for i = 1:1:TT 
	s = A.*g.^eps-c1.*g-c2.*omega-c0;
	eta = a.*m.*zeta.*g.^eps+ b.*m.*omega.^mu+f0;

	if p1>p2
		dg = (s-eta.*g).*DT;
		domega = (-eta.*g).*DT;
	
	else 
		dg = (-eta.*g).*DT;
		domega = (s-eta.*g).*DT;
	end

	
	dp1 = (p1.*(rho+eta)-eps.*zeta.*A.*g.^(eps-1)).*DT;
	dp2 = (p2.*(rho+eta)-b.*mu.*omega.^(eps-1)).*DT;

	p1 = dp1 + p1;
	p2 = dp2 + p2;
	g = dg + g;
	omega = domega + omega;

	if mod(i,TOUT) == 0
		gg(i/TOUT+1) = g;
		omegalist(i/TOUT+1) = omega;
		p1list(i/TOUT+1) = p1;
		p2list(i/TOUT+1) = p2;	
		tt(i/TOUT+1) = i*DT; 
	end
end


plot(omegalist,gg,'m');
legend('η=0','s=0','g,ω-定常状態','動的経路');

%%%%%%%%%%%%%%%%%%%PLOT%%%%%%%%%%
figure(2)
plot(tt,gg)
xlabel('時間(tt)');
ylabel('一人あたりの公共生産施設(g)')

figure(3)
plot(tt,p1list,'b'); hold on
plot(tt,p2list,'r');
legend('p1','p2');
xlabel('時間(tt)');
ylabel('潜在価格')



%figure(2);
%contour(x,y,s,[0,0],'y');
%figure(3);
%contour(x,y,h,[0,0],'r');
%}

