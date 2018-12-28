clear;
clf;

[x,y] = ndgrid(0:0.01:0.5,0:0.01:3);


%%PARAMETERS

c0 = 0.01; %その他の経常経費
c1 = 0.1; %keijo-keihiritu　gにかかる 
c2 = 0.4; % omegaにかかる経常経費

%%INTERIOR PARAMETERS

hl = 1.5; %h*l h:税率、l:拡大係数

m = 1; %移動係数.効用の差でどれくらい人が移動するか。
m1 = 1;
m2 = 1;
rho = 0.1;%主観的割引率

%%Utility Function value d*a*g^eps*zeta+b*omega^mu %%%

mu = 0.50; %μ
eps = 0.20; %ε
a = 0.625;
d = 0.1; %hl = md =0.3
b = 2.0;
zeta = 0.2;%ζ

%%%% CHUNK %%%

A = a*hl*zeta; %0.15;
B = a*d*m*zeta; %B = 0.15; %B=a*d*m*zeta
alpha = b*m; %alpha = 6.0; %α = bm 
f0 = -0.59;


%%%% plot g=omega=0 %%%
i = B.*y.^eps+b.*m.*x.^mu+f0; 
s = A.*y.^eps-c1.*y-c2.*x-c0;
h = A.*y.^eps-c1.*y-c2.*x-c0-(y+x).*(B.*y.^eps+alpha.*x.^mu+f0);

figure(1);
I = contour(x,y,i,[0,0],'b'); hold on;
S = contour(x,y,s,[0,0],'r'); hold on;
H = contour(x,y,h,[0,0],'g'); hold on;

xlabel('一人あたりの公共生活施設量(ω)');
ylabel('一人あたりの公共生産施設量(g)');


%%%%%%%% 動学経路　DIFFERENTIAL EQUATION %%%%%%%

DT = 0.001;
TT = 100000;
TOUT = 10;

%%%%% INITIAL STATE %%%%%
g = 1.2;
omega = 0.3; 
p1 = 0.3;
p2 = 0.1;  

tt(1) = 0.0;
gg(1) = g;
omegalist(1) = omega;
p1list(1) = p1;
p2list(1) = p2;

%%%% EULER METHOD %%%%

for i = 1:1:TT 
	s = A.*g.^eps-c1.*g-c2.*omega-c0;
	eta = a.*d.*m.*zeta.*g.^eps+ b.*m.*omega.^mu+f0;

	if p1*m1>p2*m2
		dg = (s*m1-eta.*g).*DT;
		domega = (-eta.*g).*DT;
	
	else 
		dg = (-eta.*g).*DT;
		domega = (s*m2-eta.*g).*DT;
	end

	
	dp1 = (p1.*(rho+eta)-eps.*a*d.*zeta.*g.^(eps-1)).*DT;
	dp2 = (p2.*(rho+eta)-b.*mu.*omega.^(mu-1)).*DT;

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


plot(omegalist,gg,'m'); hold on;

legend('η=0 (η:人口増加率)','s=0 (s:余剰財源)','g,ω-定常状態','動的経路');

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

