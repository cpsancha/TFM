load('bx_lf.mat')
load('hc_lf.mat')
load('hk_lf.mat')
load('bx_af.mat')
load('hc_af.mat')
load('hk_af.mat')
load('bottom.mat')

lx = linspace(0, AC.Hull.Lf,500);
la = linspace(AC.Hull.Lf, AC.Hull.Lf+AC.Hull.La, 500);


bx_int = interp1(bx(:,1)./4.*AC.Hull.Lf, bx(:,2).* AC.Hull.Beam,lx,'pchip');
hc_int = interp1(hc(:,1)./4.*AC.Hull.Lf, hc(:,2).* AC.Hull.Beam./0.75.*0.65,lx,'pchip');
hk_int = interp1(hk(:,1)./4.*AC.Hull.Lf, hk(:,2).* AC.Hull.Beam./0.75.*0.65,lx,'pchip');
bxa_int = interp1(bx_af(:,1)./3.5.*AC.Hull.La, bx_af(:,2).* AC.Hull.Beam,la-AC.Hull.Lf,'pchip');
hca_int = interp1(hc_af(:,1)./3.5.*AC.Hull.La, hc_af(:,2).* AC.Hull.Beam./0.75.*0.65,la-AC.Hull.Lf,'pchip');
hka_int = interp1(hk_af(:,1)./3.5.*AC.Hull.La, hk_af(:,2).* AC.Hull.Beam./0.75.*0.65,la-AC.Hull.Lf,'pchip');
p_bot  = polyfit(bottom(:,1),bottom(:,2),2);

figure()
hold on
axis equal
plot(bx(:,1)./4.*AC.Hull.Lf,bx(:,2).* AC.Hull.Beam)
plot(lx, bx_int,'r')
plot(hc(:,1)./4.*AC.Hull.Lf, hc(:,2).* AC.Hull.Beam,'g')
plot(lx, hc_int,'k')
plot(hk(:,1)./4.*AC.Hull.Lf, hk(:,2).* AC.Hull.Beam,'g')
plot(lx, hk_int,'k--')
% 
% plot(hk_af(:,1)./4.*AC.Hull.Lf, hk_af(:,2).* AC.Hull.Beam,'m')
plot(la-AC.Hull.Lf, hka_int,'m--')



stationPoints = 30;
X = [];
Y = [];
Z = [];
for i=1:(length(lx)+length(la))
    if i<=length(lx)
    x_station  = lx(i);
    hk_station = hk_int(i);
    hc_station = hc_int(i);
    bx_station = bx_int(i);
    
    y=linspace(0,1,stationPoints);
    z= polyval(p_bot,y);
    y_station = y.*bx_station./2;
    z_station = z.*(hc_station-hk_station) + hk_station;
 
    X = horzcat(X,x_station.*ones(1,2*stationPoints));
    Y = horzcat(Y,-fliplr(y_station),y_station);
    Z = horzcat(Z,fliplr(z_station),z_station);
    else
    x_station  = la(i-length(lx));
    hk_station = hka_int(i-length(lx));
    hc_station = hca_int(i-length(lx));
    bx_station = bxa_int(i-length(lx));
    
    y = linspace(0,1,stationPoints);
    y_station = y.*bx_station./2;
    z_station = y_station.*(hc_station-hk_station)./(bx_station./2) + hk_station;
    
    X = horzcat(X,x_station.*ones(1,2*stationPoints));
    Y = horzcat(Y,-fliplr(y_station),y_station);
    Z = horzcat(Z,fliplr(z_station),z_station);
    end
end

figure()
axis equal
pcshow([X(:),Y(:),Z(:)]);

malla = [X',Y',Z'];
save('malla.txt','malla','-ascii')
%%

% figure()
% hold on
% 
% plot(Data002(:,1),Data002(:,2))
% bot = polyfit(Data003(:,1),Data003(:,2),3);
% plot(linspace(0,1,100),polyval(bot,linspace(0,1,100)))
