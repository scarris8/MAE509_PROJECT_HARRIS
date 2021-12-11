%% CONSERVATIVE ROBUST METHOD
clear all; close all;

disp('System')

%Conversions
km2m = 1000;
hr2s = 3600;

%Orbital Specs
g = 1.62; %m/s2 (Lunar gravity)
M = 7.34767309e22; %kg (mass of moon)
G = 667408e-11; %m3/kg*s2
R = 1737.4; %km (radius of moon)
apo_alt = 122.4; %km (apoapsis altitude)
peri_alt = 100.9; %km (periapsis altitude)
T = 2; %hr (orbital period)
a = ((R+apo_alt)+(R+peri_alt))/2; %km (orbit semi-major axis)
Vorbit = 2*pi*(a*km2m)/(T*hr2s); %m/s (average orbital velocity)
h0 = a-R; %km (average orbital altitude above moon surface)
V0 = 0;

%Lander Specs
m_fuel = 8200; %kg
m_0 = 10334;
m_dry = m_0-m_fuel;
m_min = m_dry+m_fuel*.3;
% m = (m_0+m_min)/2;

Fmax = 45040; %N
Fmin = .1*Fmax;
F = (Fmax+Fmin)/2;
Isp = 3050; %N*s/kg

mdot = F/(Isp*g); %kg/s

%Calculate fuel necessary for orbital slowdown
xdot0 = Vorbit;
syms t
eqn = (Fmax/m_0)*t-xdot0 == 0;
S = solve(eqn,t);
S = double(S);
dt = S(S>0);
dt = min(dt);

m_loss_orbit = 1/(Isp*g)*Fmax*dt;
m_0_drop = m_0-m_loss_orbit;

dt = 1;
% t = 0:dt:1000;

%Initial Conditions
x0 = h0*km2m;
xdot0 = 0;
% x = zeros(2,length(t));

%Bounding Conditions Calculations
%Upper Bound
mmax = m_0_drop;
mmin = 4000;
tol = .001;
while mmax-mmin > 1
    m = [];
    m(1,1) = .5*(mmax+mmin);
    Fmax = 45040; %N
    Fmin = .1*Fmax;
    while Fmax-Fmin > 1
        F = .5*(Fmin+Fmax);

        x = [];
        x(1,1) = x0;
        x(2,1) = xdot0;
        x(1,2) = x0+xdot0*dt;

        
        i = 1;
        while x(1,i+1) > tol
            if m(i) > m_dry
                m(i+1) = m(i)-1/(Isp*g)*F*dt;
                x(1,i+2) = 2*x(1,i+1)-x(1,i)+(F/m(i+1)-g)*dt*dt;
            else
                m(i+1) = m(i);
                x(1,i+2) = 2*x(1,i+1)-x(1,i)+(-g)*dt*dt;
            end
            x(2,i+2) = x(1,i+2)-x(1,i+1);
            if x(2,i+2) > 0
                break
            end
            i = i+1;
        end
        if x(2,end) < -tol
            Fmin = F;
        end
        if x(2,end) > tol
            Fmax = F;
        end
    end
    if m(1,end) <= m_dry+tol*m_dry || x(2,end) >= 0
        mmin = m(1);
    else
        mmax = m(1);
    end
end

F_min_case = F;
x_min_case = x;
m_min_case = m;
t_min_case = 0:dt:dt*(length(x_min_case)-1);

%LOWER BOUND
m = [];
m(1) = m_0_drop;
Fmax = 45040; %N
Fmin = .1*Fmax;
while Fmax-Fmin > 1
    F = .5*(Fmin+Fmax);

    x = [];
    x(1,1) = x0;
    x(2,1) = xdot0;
    x(1,2) = x0+xdot0*dt;

    i = 1;
    while x(1,i+1) > tol
        if m(i) > m_dry
            m(i+1) = m(i)-1/(Isp*g)*F*dt;
            x(1,i+2) = 2*x(1,i+1)-x(1,i)+(F/m(i+1)-g)*dt*dt;
        else
            m(i+1) = m(i);
            x(1,i+2) = 2*x(1,i+1)-x(1,i)+(-g)*dt*dt;
        end
        x(2,i+2) = x(1,i+2)-x(1,i+1);
        if x(2,i+2) > 0
            break
        end
        i = i+1;
    end
    if x(2,end) < -tol
        Fmin = F;
    end
    if x(2,end) > tol
        Fmax = F;
    end
end

F_max_case = F;
x_max_case = x;
m_max_case = m;
t_max_case = 0:dt:dt*(length(x_max_case)-1);

%--------------------------------
%SYSTEM
%--------------------------------

%Nominal System
A = [0 1; 0 0];
B1 = [0; 0];
B2 = [0; 0];
C1 = [0 0];
C2 = [1 0];
D11 = 0;
D12 = 0;
D21 = 0;
D22 = 0;

n = size(A,1);

m = [];
segments = 1:length(m_max_case)/11:length(m_max_case);
segments(end+1) = length(m_max_case);
for i = 1:length(segments)-1
    m(1,i) = m_max_case(segments(i));
    m(2,i) = m_min_case(segments(i+1));
end

p = length(m);

%Plot Bounding Conditions
subplot(2,1,1);
grid on
hold on
plot(segments(1:end-1),m,'r')
for i = 1:length(segments)-1
    plot([segments(i) segments(i)],m(:,i),'b');
end
ylabel('Lander Mass (kg)');
xlabel('Time (s)');

%%
disp('LMI')

dt = 1;
%Nominal System
A = [0 1; 0 0];
B = [0; 0];
C = [0 0];
D = 0;

n = size(A,1);
p = length(m);

%Uncertain System
for i = 1:p             
    Ai(:,:,i) = [0 0; 0 0];
    Bi(:,:,i) = [0; 1/m(i)];
end

%-------------------------------
%YALMIP SECTION
%-------------------------------

yalmip('clear');
P = sdpvar(n);
Z = sdpvar(1,n,'full');
G = sdpvar(n,n,'full');
eta = 0.0001;
Const = [];
Const = [Const, P >= eta*eye(size(P))];
p_os = 10;
c=(log(p_os)/pi);
FF = zeros(p,2);
for i = 1:p
    for j = 1:2
        %Uncertain System
        Ai(:,:,i) = [0 0; 0 0];
        Bi(:,:,i) = [0; 1/m(j,i)];
        %Combine Nominal and Uncertain Systems
        AA(:,:,i) = A+Ai(:,:,i);
        BB(:,:,i) = B+Bi(:,:,i);
        CC = C;
        DD = D;
        %Discretize the System
        sysc = ss(AA(:,:,i),BB(:,:,i),CC,DD);
        sysd = c2d(sysc,dt);
        
        Ad = sysd.a;
        Bd = sysd.b;
        Cd = sysd.c;
        Dd = sysd.d;

        %Discrete Time Stabilizability
        M11 = P;
        M12 = Ad*P+Bd*Z;
        M21 = M12';
        M22 = P;
        M = [M11 M12; M21 M22];
        Const = [Const, M >= eta*eye(size(M))];

        %Robust Stability Conditions
        M11 = P;
        M12 = Ad*G-Bd*Z;
        M21 = M12';
        M22 = G+G'+P;
        M = [M11 M12; M21 M22];
        Const = [Const, M >= eta*eye(size(M))];
    end    
end

opt=sdpsettings('solver','sedumi','verbose',0);
optimize(Const, [], opt);

ZZ = value(Z);
PP = value(P);
GG = value(G);
FF = -ZZ*pinv(GG);

%%
disp('simulation')
g = 1.62;
Fmax = 45000;
m_fuel = 8200; %kg
m_0 = 10334;
m_dry = m_0-m_fuel;
x = [];
x0 = 10000;
xdot0 = 0;
x(1,1) = x0;
x(2,1) = xdot0;
x(1,2) = x0+xdot0*dt;
msim = [];
msim(1,1) = sum(m_min_case(1)+m_max_case(1))/2;

% dt = .1;
j = 1;
q = 1000;
for i = 1:dt:q-1
    Ai(:,:,i) = [0 0; 0 0];
    Bi(:,:,i) = [0; 1/msim(i)];

    AA(:,:,i) = A+Ai(:,:,i);
    BB(:,:,i) = B+Bi(:,:,i);
    CC = C;
    DD = D;
    
    sysc = ss(AA(:,:,i),BB(:,:,i),CC,DD);
    sysd = c2d(sysc,dt);

    Ad = sysd.a;
    Bd = sysd.b;
    Cd = sysd.c;
    Dd = sysd.d;
    
    u(i) = FF*x(:,i);
    Ft(i) = u(i)+msim(i)*g;
    if Ft(i) < 0
        Ft(i) = 0;
        u(i) = Ft(i)-msim(i)*g;
    end
    if Ft(i) > Fmax
        Ft(i) = Fmax;
        u(i) = Ft(i)-msim(i)*g;
    end
    if msim(i) < m_dry
        msim(i) = m_dry;
    end
%     %Correction for impacting surface
%     if x(1,i) < 0
%         x(1,i) = 0;
%         x(2,i) = 0;
%     end
%     if x(1,i) == 0
%         Ft(i) = 0;
%         u(i) = Ft(i)-msim(i)*g;
%     end
%     %End Correction
    x(:,i+1) = Ad*x(:,i)+Bd*u(i);
    msim(i+1) = msim(i)-1/(Isp*g)*Ft(i)*dt;
end
%Plot Simulation Solution
figure();
subplot(2,1,1);
hold on
yyaxis left
plot(x(1,:));
ylabel('Position (m)');
grid on;
yyaxis right
plot(x(2,:));
ylabel('Velocity (m/s)');
xlabel('Time (s)');
grid on;
title('Controlled System States')

subplot(2,1,2);
hold on
yyaxis left
plot(Ft);
ylabel('Thrust (N)');
xlabel('Time (s)');
grid on;
yyaxis right
plot(msim);
ylabel('System Mass (kg)');
title('Thrust');