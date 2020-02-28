P1=[0.000000e+00, 0.000000e+00, 0.000000e+00];
P2=[2.923852e-01, -1.844820e-01, 9.383375e-01];
P3=[-6.604095e-01, -1.567369e-01, 1.240682e+00];
P4=[-9.527947e-01, 2.774510e-02, 3.023449e-01];
P5=[8.181147e-02, 9.824442e-01, 1.676613e-01];
P6=[3.741967e-01, 7.979621e-01, 1.105999e+00];
P7=[-5.785981e-01, 8.257072e-01, 1.408344e+00];
P8=[-8.709833e-01, 1.010189e+00, 4.700061e-01];
C1=[5.697592e-01, 3.159923e-01, 1.037074e+00];
C2=[4.211637e-01, 5.624242e-01, 7.937364e-02];
R=5.000000e-01;
E0= [0,0,1];E1= [1,1,0]; E2= [-1,1,0]; E3= [-1,-1,0]; E4= [1,-1,0];
E=[E0;E1;E2;E3;E4];
PyramidZ = E0(3);
PyramidBase = [min(E(:, 1)), max(E(:, 1)), min(E(:, 2)), max(E(:,2))];
N = 200000;
I = 0;
temp=zeros(N,3);
%Cuboid
V1 = P4 - P1;
V2 = P2 - P1;
V3 = P5 - P1;
V1x = V1 .* rand(N,1);
V2y = V2 .* rand(N,1);
V3z = V3 .* rand(N,1);
Cuboid = V1x+V2y+V3z+P1;
Cuboid_Volume=norm(V1)*norm(V2)*norm(V3);
plot3(Cuboid(:,1),Cuboid(:,2),Cuboid(:,3),'.G'); hold on;

%Pyramid
Px=norm(E2-E1);
Py=norm(E4-E1);
mid_a=(E2+E4)/2;
mid_b=(E1+E3)/2;
P_height=norm(E0-(mid_b-mid_a));
n = nthroot(rand(N,1),3);
x = (-n+rand(N,1).*(n*Px));
y = (-n+rand(N,1).*(n*Py));
z = (1-n)*P_height;
Pyramid=[x y z];
Pyramid_volume=(1/3)*Px*Py*P_height;
plot3(Pyramid(:,1),Pyramid(:,2),Pyramid(:,3),'K');hold on;
%Cylinder
axialvec = C2 - C1;
axialvec = axialvec/norm(axialvec);
axialpoints = C1 + (C2 - C1).*rand(N,1);
circr = sqrt(rand(N,1))*R;
circtheta = rand(N,1)*2*pi;
circpoints =[cos(circtheta).*circr,sin(circtheta).*circr];
axnull = null(axialvec);
Cylinder = axialpoints + circpoints*axnull.';
Cylinder_Volume=pi*(R.^2)*norm(axialvec);
plot3(Cylinder(:,1),Cylinder(:,2),Cylinder(:,3),'.B'); hold on;
for i = 1 : N
    %if point exist on Pyramid
    Cube=[Cuboid(i,1),Cuboid(i,2),Cuboid(i,3)];
    S = PyramidBase.* (PyramidZ - Cube(:, 3)) / PyramidZ;
    if (Cube(:, 3) >= 0 && Cube(:, 3) <= PyramidZ)&&(Cube(:, 1) >= S(:, 1) && Cube(:, 1) <= S(:, 2)) && (Cube(:, 2) >= S(:, 3) && Cube(:, 2) <= S(:, 4))
        %if it also exist on Cylinder
        if dot((Cube-C1),(C2-C1))>=0 && dot((Cube-C2),(C2-C1))<=0
            if norm(cross((Cube-C1),(C2-C1)))/norm(C2-C1)<=R
                I=I+1; %count the intersection
                temp(I,:)=[Cuboid(i,1),Cuboid(i,2),Cuboid(i,3)];
            end
        end
    end
end
plot3(temp(1:I,1),temp(1:I,2),temp(1:I,3),'.r');hold on
hold off;
axis equal;
Result = (I/N)*Cuboid_Volume;
