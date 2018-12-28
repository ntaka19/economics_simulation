clear;
clf;

[x,y] = ndgrid(0:0.01:1,0:0.01:10);


%%PARAMETERS

c0 = 0.01; %その他の経常経費
c1 = 0.1; %　gにかかる経常経費 
c2 = 0.4; % omegaにかかる経常経費

%%INTERIOR PARAMETERS

hl = 2.4; %h*l h:税率、l:拡大係数
m = 1; %移動係数.効用の差でどれくらい人が移動するか。
m1 = 1;
m2 = 1;
rho = 0.04;%主観的割引率

%%効用関数 d*a*g^eps*zeta+b*omega^mu %%%

mu = 0.25; %μ
eps = 0.20; %ε

a = 0.3;
d = 0.9; %hl = md =0.3
b = 0.5;  %%この値によって、グラフの形状が大きく変わるのはなぜ？1.1だと接点を持たない。。
zeta = 0.5;

%%%% CHUNK(まとめたもの) %%%

A = a*hl*zeta; %0.15;
B = a*d*m*zeta; %B = 0.15; %B=a*d*m*zeta
alpha = b*m; %alpha = 6.0; %α = bm 
f0 = -0.5;


%%%% plot g=omega=0 %%%
i = B.*y.^eps+b.*m.*x.^mu+f0; 
s = A.*y.^eps-c1.*y-c2.*x-c0;
%h = A.*y.^eps-c1.*y-c2.*x-c0-(y+x).*(B.*y.^eps+alpha.*x.^mu+f0);
h = s-(y+x).*i;

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

%%%%% 初期値 %%%%%
g = 2.0;
omega = 0.3; 
p1 = 0.15;
p2 = 2.1;
u = 0;
u_sum = 0; 
du_sum_ger = 0;  

tt(1) = 0.0;
gg(1) = g;
omegalist(1) = omega;
p1list(1) = p1;
p2list(1) = p2;
du_sum_gerlist(1) = du_sum_ger*exp(-rho*0);


%%%% EULER METHOD %%%%
	for i = 1:1:TT 
		s = A.*g.^eps-c1.*g-c2.*omega-c0;
		eta = a.*d.*m.*zeta.*g.^eps+ b.*m.*omega.^mu+f0;

		u_before = u;
		u_after = a.*d.*zeta.*g.^eps+b.*omega.^mu;
		u = u_after;


%%%%潜在価格が0で打ち切り%%%%
		if p1<0 || p2 <0 
			break;
		end

		if p1*m1>=p2*m2
			dg = (s*m1-eta.*g).*DT;
			domega = (-eta.*omega).*DT;
		
		else 
			dg = (-eta.*g).*DT;
			domega = (s*m2-eta.*omega).*DT;
		end

		
		dp1 = (p1.*(rho+eta)-eps.*a*d.*zeta.*g.^(eps-1)).*DT;
		dp2 = (p2.*(rho+eta)-b.*mu.*omega.^(mu-1)).*DT;

		p1 = dp1 + p1;
		p2 = dp2 + p2;
		g = dg + g;
		omega = domega + omega;


		if p1<0 || p2 <0 
			break;
		else 
			if mod(i,TOUT) == 0
				gg(i/TOUT+1) = g;
				omegalist(i/TOUT+1) = omega;
				p1list(i/TOUT+1) = p1;
				p2list(i/TOUT+1) = p2;	
				etalist(i/TOUT+1) =eta;
				slist(i/TOUT+1) = s;
				du_sum_gerlist(i/TOUT+1) = du_sum_ger;
				ulist(i/TOUT+1) = u_before;
				tt(i/TOUT+1) = i*DT; 
			end
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
plot(tt,omegalist)
xlabel('時間(tt)');
ylabel('一人あたりの公共生活施設(ω)')

figure(4)
plot(tt,p1list,'b'); hold on
plot(tt,p2list,'r');
legend('p1','p2');
xlabel('時間(tt)');
ylabel('潜在価格')


break;
saveas(figure(1),'data4-1','png')
saveas(figure(2),'data4-2','png')
saveas(figure(3),'data4-3','png')
saveas(figure(4),'data4-4','png')


data = horzcat(tt',p1list',p2list',gg',omegalist',etalist',slist');
csvwrite('C:\Users\成章\Desktop\economics\economics\matlab_simulation\data4\data4.csv',data)


break;
fprintf(fileID,'%f,%f\n',tt',p1list)
fprintf(fileID1,'%f,%f,%f\n',tt,p1list,p2list)

fileID = fopen('pp1.txt','w')
fileID1 = fopen('pp2.txt','w')