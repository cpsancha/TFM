
run('NACA63A210_Coordinates.m')

xpos=AC.Wing1.Root_LE;
zpos=AC.Fuselage.fusHeight;

clear X Y Z x y z
X = [];
Y = [];
Z = [];
for i=1:length(AC.Wing1.eta)
    y_local = AC.Wing1.eta(i)*AC.Wing1.WingSpan/2;
    x = xpos+NACA63A210(:,1).*AC.Wing1.c(i)+ (AC.Wing1.RootChord-AC.Wing1.c(i))/2;
    z = zpos+NACA63A210(:,2).*AC.Wing1.c(i);
    
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
Wing1Mesh = [X',Y',Z'];
save('Wing1Mesh.txt','Wing1Mesh','-ascii')
