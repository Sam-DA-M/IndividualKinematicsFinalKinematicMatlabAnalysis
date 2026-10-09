%Kinematics of Mechanisms: Dynamic Analysis of Wiper Mechanism.
%Grounded Joints:AD

%Joint Coordinates (in meters)

A = [-.938 0.476 0];
AAdjusted = A*.081753
A = AAdjusted
B = [-1.49 .39 0];
BAdjusted = B*.081753;
B = BAdjusted
C = [.517 -.119 0];
CAdjusted = C*.081753;
C=CAdjusted
D = [.87 .697 0];
DAdjusted = D*.081753;
D = DAdjusted
E = [3.79 1 0];
EAdjusted = E*.081753;
E = EAdjusted

%Link Lengths

% Define link lengths based on the coordinates
lengthAB = norm(A-B)
lengthBC = norm(B-C)
lengthCD = norm(D-C)
lengthDE = norm(E-D)
lengthEC = norm(E-C)

%Angle EC:

vectorDC = C - D
vectorDE = E - D
ThetaEC = acos(dot(vectorDC,vectorDE)/(lengthDE*lengthCD));

%Other dimensions
width = .02;
thickness = .012;
jointdiameter = .06;
density = 7870;
%Other material properties
% E = 2.05e+11 B/m^2
% poison's = .29
%Tensile strength = 385000000 Pa
%Yield strength = 325000000

% Define link masses
massAB = .10115;
massBC = .34944;
massDE = .618;

%Define link weights
WAB = [0 -massAB*9.81 0]
WBC = [0 -massBC*9.81 0]
WDE = [0 -massDE*9.81 0]


%Center of mass of each link
S1 = (A+B)/2;
S2 = (B+C)/2;
S3estimate = [1.8 .5 0];
S3 = S3estimate*.081753;

lengthS3D = norm(S3-D)
lengthS3E = norm(S3-E)

%%




%Newton's II law
%Angular velocity and angular acceleration?
%Calculating instantaneous angular velocity under assumptions.
omegaAB_RPM = 75 %RPM
omegaAB = omegaAB_RPM*(2*pi/60) %rad/s

%%

%Angular Acceleration
%Calculating instantaneous angular acceleration under assumptions.
%alphaGF = d_omega/d_t = (d_theta/d_t) * (d_omega/d_theta)
alphaAB = 0

vector_alphaAB = [0 0 alphaAB];


%%

%Newton's Second Law
%TBD!!!!
%complex mass moment of inertia
J_AB = 32698.21/(1000^3);
J_ED = 5022784.92/(1000^3);
J_CB = 1005339.39/(1000^3);
%%


%Positional Analysis

%Circle Intersection Technique

%Joint coordinates have been defined
%length of links also defined.

%Compute initial angle of input link AB.

initial_theta = atan2(B(2)-A(2),B(1)-A(1));

if (initial_theta<0)
   inputAngle = 2*pi + initial_theta;
else
    inputAngle = initial_theta;
end

for theta=0:1:360
% new position of joint A
theta;

B_new = [A(1)+lengthAB*cos(inputAngle+deg2rad(theta)), A(2)+lengthAB*sin(inputAngle+deg2rad(theta)), 0];

[Cx,Cy]= circcirc(B_new(1),B_new(2),lengthBC,D(1),D(2),lengthCD);

%Check if NaN

circIntersect_x = any(isnan(vpa(Cx)));
circIntersect_y = any(isnan(vpa(Cy)));

if circIntersect_x == 0 && circIntersect_y == 0
C_1 = [Cx(1) Cy(1) 0];
C_2 = [Cx(2) Cy(2) 0];

%distance to determine whether C_1 or C_2 is correct

dist1 = norm(C_1-C);
dist2 = norm(C_2-C);

if(dist1<dist2)
    C_new = vpa(C_1);

else

    C_new = vpa(C_2);
end

%New Position of Joint E using C_new

[Ex,Ey]= circcirc(C_new(1),C_new(2),lengthEC,D(1),D(2),lengthDE);

%Check if there is NaN
circIntersect_x = any(isnan(vpa(Ex)));
circIntersect_y = any(isnan(vpa(Ey)));

if circIntersect_x == 0 && circIntersect_y == 0
    E_1 = [Ex(1) Ey(1) 0];
    E_2 = [Ex(2) Ey(2) 0];


    dist1 = norm(E_1-E);
    dist2 = norm(E_2-E);


if(dist1<dist2)
    E_new = vpa(E_1);

else

    E_new = vpa(E_2);
end


%New Position of S3 using E_new and D
[S3x,S3y]= circcirc(E_new(1),E_new(2),lengthS3E,D(1),D(2),lengthS3D);

%Check if there is NaN
circIntersect_x = any(isnan(vpa(S3x)));
circIntersect_y = any(isnan(vpa(S3y)));

if circIntersect_x == 0 && circIntersect_y == 0
    S3_1 = [S3x(1) S3y(1) 0];
    S3_2 = [S3x(2) S3y(2) 0];


    dist1 = norm(S3_1-S3);
    dist2 = norm(S3_2-S3);

if(dist1<dist2)
    S3_new = vpa(S3_1);

else

    S3_new = vpa(S3_2);
end



%Store values for plotting
new_B_x((theta+1),:) = B_new(1);
new_B_y((theta+1),:) = B_new(2);
new_C_x((theta+1),:) = C_new(1);
new_C_y((theta+1),:) = C_new(2);
new_E_x((theta+1),:) = E_new(1);
new_E_y((theta+1),:) = E_new(2);

B = B_new;
C = C_new;
E = E_new;
S3=S3_new;

else
    fprintf('New position of S3 cannot be determined at angle: %d degree', theta)
break
end

else
    fprintf('New position of E cannot be determined at angle: %d degree', theta)
break
end

else
    fprintf('New position of C cannot be determined at angle: %d degree', theta)
break
end



%Implicit Center of mass of links
S1 = (A+B)/2;
S2 = (B+C)/2;
S3_PAtest((theta+1),:)=S3_new;

%Calculate forces and torque based on static equilibrium conditions
syms FAx FAy FBx FBy FCx FCy FDx FDy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
InputTorque = [0 0 Tin];


%Forces
%ForceAB
eqn1 = ForceA + ForceB + WAB == 0;
%ForceBC
eqn2 = -ForceB + ForceC + WBC == 0;
%ForceCDE
eqn3 = -ForceC + ForceD + WDE == 0;

%Moments
%MomentAB
eqn6 = InputTorque + cross(A-S1,ForceA) + cross(B-S1,ForceB) == 0;
%MomentBC
eqn7 = cross(C-S2, ForceC) + cross(B-S2,-ForceB) == 0;
%MomentCDE
eqn8 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) == 0;

eqnmatrix = [eqn1, eqn2, eqn3, eqn3, eqn6, eqn6, eqn7, eqn8];
% Solve the system of equations for the unknown forces and torque
staticsolution = solve(eqnmatrix, [FAx, FAy, FBx, FBy, FCx, FCy, FDx, FDy, Tin]);

%storing values
x_theta((theta+1),1)=theta;
ForceA_PA((theta+1),:) = sqrt((double(staticsolution.FAx))^2 + (double(staticsolution.FAy))^2);
ForceB_PA((theta+1),:) = sqrt((double(staticsolution.FBx))^2 + (double(staticsolution.FBy))^2);
ForceC_PA((theta+1),:) = sqrt((double(staticsolution.FCx))^2 + (double(staticsolution.FCy))^2);
ForceD_PA((theta+1),:) = sqrt((double(staticsolution.FDx))^2 + (double(staticsolution.FDy))^2);

vectorForceA_PA((theta+1),:) = [double(staticsolution.FAx) double(staticsolution.FAy), 0];
vectorForceB_PA((theta+1),:) = [double(staticsolution.FBx) double(staticsolution.FBy), 0];
vectorForceC_PA((theta+1),:) = [double(staticsolution.FCx) double(staticsolution.FCy), 0];
vectorForceD_PA((theta+1),:) = [double(staticsolution.FDx) double(staticsolution.FDy) 0];
vectorTorquein_PA((theta+1),:) = [0 0 double(staticsolution.Tin)];



%Newton's II law
%Angular velocity and angular acceleration?


%loop ADCBA
syms wCB wED

vector_omegaAB = [0 0 omegaAB];
omegaCB = [0 0 wCB];
omegaED = [0 0 wED];

eqn12 = cross(omegaED, D-C) + cross(omegaCB, C-B) + cross(vector_omegaAB, B-A) == 0;
loop2Solution = solve(eqn12,[wCB,wED]);

angularVelocity_CB = double(loop2Solution.wCB);
angularVelocity_ED = double(loop2Solution.wED);

vector_omegaCB = [0 0 angularVelocity_CB];
vector_omegaED = [0 0 angularVelocity_ED];

vector_omegaCB_PA((theta+1),:) = [0 0 angularVelocity_CB];
vector_omegaAB_PA((theta+1),:) = vector_omegaAB;
vector_omegaED_PA((theta+1),:) = [0 0 angularVelocity_ED];


%Angular Acceleration
%Calculating instantaneous angular acceleration under assumptions.
%alphaGF = d_omega/d_t = (d_theta/d_t) * (d_omega/d_theta)
alphaAB = 0;

vector_alphaAB = [0 0 alphaAB];
vector_alphaAB_PA((theta+1),:) = [0 0 alphaAB];


%loop ADCBA

syms aCB aED

alphaCB = [0 0 aCB];
alphaED = [0 0 aED];

a_ED = cross(alphaED,D-C) + cross(vector_omegaED,cross(vector_omegaED,D-C));
a_CB = cross(alphaCB,C-B) + cross(vector_omegaCB,cross(vector_omegaCB,C-B));
a_BA = cross(vector_alphaAB,B-A) + cross(vector_omegaAB,cross(vector_omegaAB,B-A));
eqn14 = a_ED + a_CB + a_BA == 0;

loop2angularaccels = solve(eqn14,[aCB,aED]);

vector_alphaCB = [0 0 double(loop2angularaccels.aCB)];
vector_alphaED = [0 0 double(loop2angularaccels.aED)];
vector_alphaCB_PA((theta+1),:) = [0 0 double(loop2angularaccels.aCB)];
vector_alphaED_PA((theta+1),:) = [0 0 double(loop2angularaccels.aED)];

%%

%Velocity at Joints with respect to ground

%Velocity at B
vectorvB_A = cross(vector_omegaAB,B-A);

vectorvB_A_PA((theta+1),:) = cross(vector_omegaAB,B-A);

%Velocity at C
vectorvC_D = cross(vector_omegaED,C-D);
vectorvC_D_PA((theta+1),:) = cross(vector_omegaED,C-D);

%Velocity at E
vectorvE_D = cross(vector_omegaED, E-D);
vectorvE_D_PA((theta+1),:) = cross(vector_omegaED, E-D);
%%

%Accelerations at Joints with respect to ground

%Accel B
vectoraB_A = cross(vector_alphaAB,B-A) + cross(vector_omegaAB,cross(vector_omegaAB,B-A));
vectoraB_A_PA((theta+1),:) = cross(vector_alphaAB,B-A) + cross(vector_omegaAB,cross(vector_omegaAB,B-A));

%AccelC
vectoraC_D = cross(vector_alphaED,C-D) + cross(vector_omegaED,cross(vector_omegaED,C-D));
vectoraC_D_PA((theta+1),:) = cross(vector_alphaED,C-D) + cross(vector_omegaED,cross(vector_omegaED,C-D));

%AccelE
vectoraE_D = cross(vector_alphaED,E-D) + cross(vector_omegaED,cross(vector_omegaED,E-D));
vectoraE_D_PA((theta+1),:) = cross(vector_alphaED,E-D) + cross(vector_omegaED,cross(vector_omegaED,E-D));

%%

%Velocities and accelerations for each cm


%V_S1/A = cross(omega_AB,S1-A)

v_S1_A = cross(vector_omegaAB, S1-A);
v_S1_A_PA((theta+1),:) = cross(vector_omegaAB, S1-A);

%v_S2_A = V_S2_B + V B_A

v_S2_B = cross(vector_omegaCB,S2-B);
v_S2_A = v_S2_B + vectorvB_A;
v_S2_A_PA((theta+1),:) = v_S2_B + vectorvB_A;

%v_S3_D
v_S3_D = cross(vector_omegaED,S3-D);
v_S3_D_PA((theta+1),:) = cross(vector_omegaED,S3-D);

%a_S1_A
a_S1_A = cross(vector_alphaAB,S1-A) + cross(vector_omegaAB,cross(vector_omegaAB,S1-A));
a_S1_A_PA((theta+1),:) = cross(vector_alphaAB,S1-A) + cross(vector_omegaAB,cross(vector_omegaAB,S1-A));

%a_S2_A= vectoraS2_B+vectoraBA
a_S2_A = vectoraB_A + cross(vector_alphaCB,S2-B) + cross(vector_omegaCB,cross(vector_omegaCB,S2-B));
a_S2_A_PA((theta+1),:) = vectoraB_A + cross(vector_alphaCB,S2-B) + cross(vector_omegaCB,cross(vector_omegaCB,S2-B));

%a_S3_D
a_S3_D = cross(vector_alphaED,S3-D) + cross(vector_omegaED,cross(vector_omegaED,S3-D));
a_S3_D_PA((theta+1),:) = cross(vector_alphaED,S3-D) + cross(vector_omegaED,cross(vector_omegaED,S3-D));

%Newton's Second Law

syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NTin
%Define forces

NForceA = [NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD= [NFDx NFDy 0];
Torqueinput = [0 0 NTin];

%Equations for link AB
%Sum of forces
eqn15 = NForceA + NForceB + WAB == massAB*a_S1_A;
%Sum of moments
eqn16 = cross(A-S1,NForceA) + cross(B-S1,NForceB) + Torqueinput == J_AB*vector_alphaAB;

%Equations for link BC
%Sum of forces
eqn17 = -NForceB + NForceC + WBC == massBC*a_S2_A;
%Sum of moments
eqn18 = cross(B-S2,-NForceB) + cross(C-S2,NForceC) == J_CB*vector_alphaCB;

%Equations for link CDE
%Sum of forces
eqn19 = -NForceC + NForceD + WDE == massDE*a_S3_D;
%Sum of moments
eqn20 = cross(C-S3,-NForceC) + cross(D-S3,NForceD) == J_ED*vector_alphaED;

%Solving equations
NeqnMatrix = [eqn15, eqn16, eqn17, eqn18, eqn19, eqn20];
DynamicSolution = solve(NeqnMatrix, [NFAx, NFAy, NFBx, NFBy, NFCx, NFCy, NFDx, NFDy, NTin]);

%Extract forces from the dynamic solution
NForce_Ax = double(DynamicSolution.NFAx);
NForce_Ay = double(DynamicSolution.NFAy);
NForce_A = sqrt((NForce_Ay)^2+(NForce_Ax)^2);
NForce_Bx = double(DynamicSolution.NFBx);
NForce_By = double(DynamicSolution.NFBy);
NForce_B = sqrt((NForce_By)^2+(NForce_Bx)^2);
NForce_Cx = double(DynamicSolution.NFCx);
NForce_Cy = double(DynamicSolution.NFCy);
NForce_C = sqrt((NForce_Cy)^2+(NForce_Cx)^2);
NForce_Dx = double(DynamicSolution.NFDx);
NForce_Dy = double(DynamicSolution.NFDy);
NForce_D = sqrt((NForce_Dy)^2+(NForce_Dx)^2);
NTorque = double(DynamicSolution.NTin);

VectorNTorque = [0, 0, NTorque];

VectorNForceA_PA((theta+1),:) = NForce_A;

VectorNForceB_PA((theta+1),:) = NForce_B;

VectorNForceC_PA((theta+1),:) = NForce_C;

VectorNForceD_PA((theta+1),:) = NForce_D;

VectorNTorque_PA((theta+1),:) = [0, 0, NTorque];

NForceA_PA((theta+1),:) = [NForce_Ax NForce_Ay 0];
NForceB_PA((theta+1),:) = [NForce_Bx NForce_By 0];
NForceC_PA((theta+1),:) = [NForce_Cx NForce_Cy 0];
NForceD_PA((theta+1),:) = [NForce_Dx NForce_Dy 0];
end


%%
%Plots
plot(x_theta,ForceA_PA(:,1),"r-"); hold on

plot(x_theta,ForceB_PA(:,1),"y-")

plot(x_theta,ForceC_PA(:,1),"g-")

plot(x_theta,ForceD_PA(:,1),"c-")

plot(x_theta,vectorTorquein_PA(:,3),"r-.")
legend("FA","FB", "FC","FD","Tin")
title("Static Forces vs Angle")
xlabel('theta(degrees)')
ylabel('Force(N)')
hold off

%%



plot(x_theta,vector_omegaAB_PA(:,3),"r-"); hold on
plot(x_theta,vector_omegaCB_PA(:,3),"y-")
plot(x_theta,vector_omegaED_PA(:,3),"g-")
legend("omegaAB","omegaCB","omegaED")
title("Link Angular Velocity vs Angle")
xlabel('theta(degrees)')
ylabel('Angular Velocity (rad/s)')
hold off

figure(3);
plot(x_theta,vector_alphaAB_PA(:,3),"r-"); hold on
plot(x_theta,vector_alphaCB_PA(:,3),"y-")
plot(x_theta,vector_alphaED_PA(:,3),"g-")
legend("alphaAB","alphaCB","alphaED")
xlabel('theta(degrees)')
ylabel('Link Angular Acceleration (rad/s^2)')
title('Angular Acceleration vs Angle')

hold off


plot(x_theta,vectorvB_A_PA(:,1),"r-"); hold on

plot(x_theta,vectorvB_A_PA(:,2),"r--")

plot(x_theta,vectorvC_D_PA(:,1),"y-")

plot(x_theta,vectorvC_D_PA(:,2),"y--")

plot(x_theta,vectorvE_D_PA(:,1),"g-")

plot(x_theta,vectorvE_D_PA(:,2),"g--")
legend("vBx", "vBy","vCx", "vCy","vEx", "vEy")
title("Velocity at Joints with respect to ground vs Angle")
xlabel('theta(degrees)')
ylabel('Velocity (m/s)')
hold off

figure(5);
plot(x_theta,vectoraB_A_PA(:,1),"r-"); hold on

plot(x_theta,vectoraB_A_PA(:,2),"r--")

plot(x_theta,vectoraC_D_PA(:,1),"y-")

plot(x_theta,vectoraC_D_PA(:,2),"y--")

plot(x_theta,vectoraE_D_PA(:,1),"g-")

plot(x_theta,vectoraE_D_PA(:,2),"g--")
legend("aBx", "aBy","aCx", "aCy","aEx", "aEy")
title("Acceleration at Joints with respect to ground vs Angle")
xlabel('theta(degrees)')
ylabel('Acceleration (m/s^2)')
hold off


plot(x_theta,v_S1_A_PA(:,1),"r-"); hold on

plot(x_theta,v_S1_A_PA(:,2),"r--")

plot(x_theta,v_S2_A_PA(:,1),"y-")

plot(x_theta,v_S2_A_PA(:,2),"y--")

plot(x_theta,v_S3_D_PA(:,1),"g-")

plot(x_theta,v_S3_D_PA(:,2),"g--")
legend("vS1x", "vS1y","vS2x", "vS2y","vS3x", "vS3y")
title("Velocity at CM's with respect to ground vs Angle")
xlabel('theta(degrees)')
ylabel('Velocity (m/s)')
hold off

plot(x_theta,a_S1_A_PA(:,1),"r-"); hold on

plot(x_theta,a_S1_A_PA(:,2),"r--")

plot(x_theta,a_S2_A_PA(:,1),"y-")

plot(x_theta,a_S2_A_PA(:,2),"y--")

plot(x_theta,a_S3_D_PA(:,1),"g-")

plot(x_theta,a_S3_D_PA(:,2),"g--")
legend("aS1x", "aS1y","aS2x", "aS2y","aS3x", "aS3y")
title("Acceleration at CM's with respect to ground vs Angle")
xlabel('theta(degrees)')
ylabel('acceleration (m/s^2)')
hold off

plot(x_theta,VectorNForceA_PA(:,1),"r-"); hold on

plot(x_theta,VectorNForceB_PA(:,1),"y-")

plot(x_theta,VectorNForceC_PA(:,1),"g-")

plot(x_theta,VectorNForceD_PA(:,1),"c-")

plot(x_theta,VectorNTorque_PA(:,3),"r:")
legend("FA","FB","FC","FD","Tin")
title("Newton's Second Law vs Angle")
xlabel('theta(degrees)')
ylabel('Force (N)')
hold off

plot(new_B_x(:,1),new_B_y(:,1),"y-"); hold on

plot(new_C_x(:,1),new_C_y(:,1),"g-")

plot(new_E_x(:,1),new_E_y(:,1),"m-")

legend("B","C","E")
title("y Position of Joints vs x Position of Joints")
xlabel('x Position (m)')
ylabel('y Position (m)')
hold off



%%
%This section reserved for table-making and data evaluation.
forcesmatrix(1,1:3)=vectorForceA_PA(1,1:3);
forcesmatrix(2,1:3)=vectorForceB_PA(1,1:3);
forcesmatrix(3,1:3)=vectorForceC_PA(1,1:3);
forcesmatrix(4,1:3)=vectorForceD_PA(1,1:3);
forcesmatrix(5,1:3)=vectorTorquein_PA(1,1:3);

forces_table=array2table(forcesmatrix, 'VariableNames',{'x','y','z'},'RowNames', {'FA(N)','FB(N)','FC(N)','FD(N)','Tin(N)'})



%%
%These sections reserved for table-making and data evaluation.
Nmatrix(1,1:3)=NForceA_PA(1,1:3);
Nmatrix(2,1:3)=NForceB_PA(1,1:3);
Nmatrix(3,1:3)=NForceC_PA(1,1:3);
Nmatrix(4,1:3)=NForceD_PA(1,1:3);
Nmatrix(5,1:3)=VectorNTorque_PA(1,1:3);

Nmatrix_table=array2table(Nmatrix, 'VariableNames',{'x','y','z'},'RowNames', {'FA(N)','FB(N)','FC(N)','FD(N)','Tin(N*m)'})

matrix_omega(1,1:3)=vector_omegaAB_PA(1,1:3);
matrix_omega(2,1:3)=vector_omegaCB_PA(1,1:3);
matrix_omega(3,1:3)=vector_omegaED_PA(1,1:3);

omegamatrix_table=array2table(matrix_omega, 'VariableNames',{'x','y','z'},'RowNames', {'wAB(rad/s)','wCB(rad/s)','wED(rad/s)'})


matrix_alpha(1,1:3)=vector_alphaAB_PA(1,1:3);
matrix_alpha(2,1:3)=vector_alphaCB_PA(1,1:3);
matrix_alpha(3,1:3)=vector_alphaED_PA(1,1:3);


alpha_table=array2table(matrix_alpha, 'VariableNames',{'x','y','z'},'RowNames', {'alphaAB(rad/s^2)','alphaCB(rad/s^2)','alphaED(rad/s^2)'})


matrix_v(1,1:3)=vectorvB_A_PA(1,1:3);
matrix_v(2,1:3)=vectorvC_D_PA(1,1:3);
matrix_v(3,1:3)=vectorvE_D_PA(1,1:3);
matrix_v=round(matrix_v,4);
vmatrix_table=array2table(matrix_v, 'VariableNames',{'x','y','z'},'RowNames', {'vBA(m/s)','vCD(m/s)','vED(m/s)'})

matrix_a(1,1:3)=vectoraB_A_PA(1,1:3);
matrix_a(2,1:3)=vectoraC_D_PA(1,1:3);
matrix_a(3,1:3)=vectoraE_D_PA(1,1:3);
matrix_a=round(matrix_a,4);
amatrix_table=array2table(matrix_a, 'VariableNames',{'x','y','z'},'RowNames', {'aBA(m/s^2)','aCD(m/s^2)','aED(m/s^2)'})


matrix_as(1,1:3)=a_S1_A_PA(1,1:3);
matrix_as(2,1:3)=a_S2_A_PA(1,1:3);
matrix_as(3,1:3)=a_S3_D_PA(1,1:3);
matrix_as=round(matrix_as,4);
asmatrix_table=array2table(matrix_as, 'VariableNames',{'x','y','z'},'RowNames', {'aS1A(m/s^2)','aS2A(m/s^2)','aS3D(m/s^2)'})