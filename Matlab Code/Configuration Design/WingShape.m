
run('NACA63A210_Coordinates.m')

xpos=AC.Wing1.Root_LE;
zpos=AC.Fuselage.fusHeight;

clear X Y Z x y z
X = [];
Y = [];
Z = [];
for i=1:length(AC.Wing1.eta)
    y_local = AC.Wing1.eta(i)*AC.Wing1.WingSpan/2;
    x = xpos+NACA63A210e(:,1).*AC.Wing1.c(i)+ (AC.Wing1.RootChord-AC.Wing1.c(i))/2;
    z = zpos+NACA63A210e(:,2).*AC.Wing1.c(i);
    
    X = horzcat(X,x');
    Y = horzcat(Y,y_local.*ones(1,length(x)));
    Z = horzcat(Z,z');
   
end

X = horzcat(X,X);
Y = horzcat(Y,-Y);
Z = horzcat(Z,Z);
 
figure()
axis equal
pcshow([X(:),Y(:),Z(:)]);
% 
Wing1Meshe = [X',Y',Z'];
save('Wing1Meshe.txt','Wing1Meshe','-ascii')

%Intrados
xpos=AC.Wing1.Root_LE;
zpos=AC.Fuselage.fusHeight;

clear X Y Z x y z
X = [];
Y = [];
Z = [];
for i=1:length(AC.Wing1.eta)
    y_local = AC.Wing1.eta(i)*AC.Wing1.WingSpan/2;
    x = xpos+NACA63A210i(:,1).*AC.Wing1.c(i)+ (AC.Wing1.RootChord-AC.Wing1.c(i))/2;
    z = zpos+NACA63A210i(:,2).*AC.Wing1.c(i);
    
    X = horzcat(X,x');
    Y = horzcat(Y,y_local.*ones(1,length(x)));
    Z = horzcat(Z,z');
   
end

X = horzcat(X,X);
Y = horzcat(Y,-Y);
Z = horzcat(Z,Z);
 
figure()
axis equal
pcshow([X(:),Y(:),Z(:)]);
% 
Wing1Meshi = [X',Y',Z'];
save('Wing1Meshi.txt','Wing1Meshi','-ascii')

%% Ala 2
xpos=AC.Wing2.Root_LE;
zpos=AC.Fuselage.fusHeight;

clear X Y Z x y z
X = [];
Y = [];
Z = [];
for i=1:length(AC.Wing2.eta)
    y_local = AC.Wing2.eta(i)*AC.Wing2.WingSpan/2;
    x = xpos+NACA63A210e(:,1).*AC.Wing2.c(i)+ (AC.Wing2.RootChord-AC.Wing2.c(i))/2;
    z = zpos+NACA63A210e(:,2).*AC.Wing2.c(i);
    
    X = horzcat(X,x');
    Y = horzcat(Y,y_local.*ones(1,length(x)));
    Z = horzcat(Z,z');
   
end

X = horzcat(X,X);
Y = horzcat(Y,-Y);
Z = horzcat(Z,Z);
 
figure()
axis equal
pcshow([X(:),Y(:),Z(:)]);
% 
Wing2Meshe = [X',Y',Z'];
save('Wing2Meshe.txt','Wing2Meshe','-ascii')

%Intrados
xpos=AC.Wing2.Root_LE;
zpos=AC.Fuselage.fusHeight;

clear X Y Z x y z
X = [];
Y = [];
Z = [];
for i=1:length(AC.Wing2.eta)
    y_local = AC.Wing2.eta(i)*AC.Wing2.WingSpan/2;
    x = xpos+NACA63A210i(:,1).*AC.Wing2.c(i)+ (AC.Wing2.RootChord-AC.Wing2.c(i))/2;
    z = zpos+NACA63A210i(:,2).*AC.Wing2.c(i);
    
    X = horzcat(X,x');
    Y = horzcat(Y,y_local.*ones(1,length(x)));
    Z = horzcat(Z,z');
   
end

X = horzcat(X,X);
Y = horzcat(Y,-Y);
Z = horzcat(Z,Z);
 
figure()
axis equal
pcshow([X(:),Y(:),Z(:)]);
% 
Wing2Meshi = [X',Y',Z'];
save('Wing2Meshi.txt','Wing2Meshi','-ascii')