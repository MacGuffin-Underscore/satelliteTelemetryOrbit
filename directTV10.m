%% directTV10.m
clc; clearvars;
% Worry about getting MATLAB to read a TLE file later.
% DIRECTV10 = [1 31862U 07032A   20035.08310355 -.00000094  00000-0  00000+0 0  9990;
%              2 31862   0.0090 274.4608 0000490  66.6555  79.7515  1.00267821 46098];

% ~~~~~~~~~~~~~~~~~~~~~~~~
% IMPORTANT
% 1 = Earth
% 2 = directTV10
% 3 = Moon
% 4 = Sun
% ~~~~~~~~~~~~~~~~~~~~~~~~

% Get Values from TLE data
Year = 2020;                            % Julian year
Day = 035.08310355;                     % Day as a fraction
nBar = 1.00267821;                      % Mean motion (rev/day)
nDot = 2 * -.00000094^2;                % Mean motion derivative (rev/day^2)
nDoubDot = 6 * .00000E-0^3;             % Mean motion 2nd derivative (rev/day^3)
BStar = .00000E+0;                      % BStar (/ER)
e = .0000490;                           % Eccentricity
Msat = 079.7515;                        % Mean Anomoly (deg)
i = 0.0090;                            % Inclination (deg)
RAAN = 274.4608;                        % Right ascension of the ascending node (deg)
AoP = 066.6555;                         % Argument of perigee

%Other Constants
G = 6.67259e-20;                        % Gravitational constant
JD = juliandate([Year, 0, Day]);        % Julian Date (day)
DT = datetime(JD,'convertfrom','juliandate');
rEarth = 6378.135;                      % Radius of earth (km)
rho0 = 2.461e-05;                       % Atmospheric density at perigee (kg/m^2/ER)
mu = [398600;...                        % Standard gravitational parameter of earth (km^3/s^2) 
        0;...                           % PLACEHOLDER
        4904.9;...                      % Standard gravitational parameter of the moon (km^3/s^2) 
        1.3271E8];                      % Standard gravitational parameter of the sun (km^3/s^2) 
w = 7.292115146706979E-5;               % Earth rotation rate
eEarth = 0.081819221456;                % Eccentricity of earths orbit

% Calculate Additional Constants
BC = (rEarth*rho0)/(2*BStar);           % Balistic Coefficient (kg/m^2)
% BC = (1/(12.741621*BStar))            % Can also be written like this.
T = nBar^(-1)*(86400);                  % Period (s)
aBar = ((T/(2*pi))^2*mu(1))^(1/3);      % Semi-major axis (km)
h = sqrt(mu(1)*aBar*(1-e)^2);
v = Msat + (2*e-0.25*e^3)*sind(Msat)... % True Anomoly (deg)
        + 1.25*e^2*sind(2*Msat)...
        + (13/12)*e^3*sind(3*Msat); 

%% Convert Classical Orbital Elements to r and V vectors at the Epoch in the ECI Frame
% (3-1-3) Euler Angles
C = C3(AoP)*C1(i)*C3(RAAN);             % DCM: ECI to orbital frame
CT = transpose(C);                      % DCM: Orbital to ECI frame

% Initial Position Vector
r12n = h^2/(mu(1)*(1+cosd(v)));         % Magnitude of r at M (km)
r12Orb = [r12n*cosd(v);...              % r vector in orbital frame (km)
        r12n*sind(v);...
        0];                            
r12 = CT*r12Orb;                        % r vector in ECI frame (km)

% Velocity Vector
r12DotRadi = (mu(1)/h)*e*sind(v);       % Magnitude of V at M  in radial direction(km/s)
r12DotPerp = (mu(1)/h)*(1+e*cosd(v));   % Magnitude of V at M  in perpendicular direction(km/s)
gam = atand(r12DotRadi/r12DotPerp);     % Flight path angle (deg)
r12DotVec = [r12DotPerp;                % V vector in orbital frame (km/s)
           r12DotRadi; 
           0]; 
r12Dotn = norm(r12DotVec);              % Magnitude of V at v (km/s)
r12DotOrbX = r12DotRadi*cosd(v) - r12DotPerp*sind(v);
r12DotOrbY = r12DotRadi*sind(v) + r12DotPerp*cosd(v);
r12DotOrb = [r12DotOrbX;...             % V vector in orbital frame (km/s)
            r12DotOrbY;...
            0];   
r12Dot = CT*r12DotOrb;                  % V vector in ECI frame (km/s)

%% Propogate the Orbit for 30 days
X0 = [r12;r12Dot];
n = 1001;                                 % Number of steps
tSpan = linspace(JD,(JD+30),n).*(86400);  % 30 days with 1000 steps (sec)

options = odeset('RelTol', 5e-10, 'AbsTol', 5e-10);
[t, X] = ode45(@orbProp,tSpan,X0,options,mu);
rSat = [X(:,1), X(:,2), X(:,3)];

%% Change of Frame from Orbital Frame to ECEF Frame
% This goes back to using the original position and velocity vectors
CEcef = dcmeci2ecef('IAU-2000/2006',[2020 2 04 1 59 40]);
r12Ecef = CEcef*r12;
for i = 1:n         % Find r in ECEF frame over propagation
    rSatEcef(i,:) = CEcef*transpose(rSat(i,:));
end

%% Find the Longitude and Latitude
lla0 = ecef2lla([transpose(r12Ecef)], 'WGS84');     % Lat and long at initial time
latDms = degrees2dms(lla0(1));                      % Lat into DMS notation
longDms = degrees2dms(lla0(2));                     % Long into DMS notation

for i = 1:n         % Find lat and long over propagation
    lla(i,:) = ecef2lla([rSatEcef(i,:)], 'WGS84');
end
figure(3)
hold on
title('Lat vs Long')
plot(lla(:,1),lla(:,2))
hold off
%% OUTPUT
% ~~~~ Question 1 ~~~~
disp(['Q1) The Date-Time for the epoch is ',char(DT),'.']);
disp('~~~~~~~~');

% ~~~~ Question 2 ~~~~
disp(['Q2) The Julian-Date for the epoch is ',num2str(JD),'.']);
disp('~~~~~~~~');

% ~~~~ Question 3 ~~~~
disp(['Q3) The 6 classical orbital elements are:']);
        disp(['     RAAN = ',num2str(RAAN),'deg']);
        disp(['     AoP = ',num2str(AoP),'deg']);
        disp(['     i = ',num2str(i),'deg']);
        disp(['     v = ',num2str(v),'deg']);
        disp(['     e = ',num2str(e)]);
        disp(['     a = ',num2str(aBar),'km']);
disp('~~~~~~~~');
        
% ~~~~ Question 4 ~~~~
disp(['Q4) The position vector in the ECI frame is: r12 = [',num2str(transpose(r12)),']km']);
disp(['    and the velocity vector in the ECI frame is: rDot12 = [',num2str(transpose(r12Dot)),']km/s']);
disp('~~~~~~~~');

% ~~~~ Question 5 ~~~~
figure(1)
hold on
title('Orbit of directTV10 in the ECI frame')
plot3(0,0,0,'.','MarkerSize',20,'MarkerFaceColor','b')
plot3(rSat(:,1),rSat(:,2),rSat(:,3))
grid on
pbaspect([1 1 1])
daspect([1 1 1])
view(-45,10)
hold off

figure(2)
hold on
title('Orbit of directTV10 in the ECI frame (no fixed frame)')
plot3(0,0,0,'.','MarkerSize',20,'MarkerFaceColor','b')
plot3(rSat(:,1),rSat(:,2),rSat(:,3))
grid on
view(-45,10)
hold off

% ~~~~ Question 6 ~~~~
disp(['Q6) The position vector in the ECEF frame is: r12Ecef = [',num2str(transpose(r12Ecef)),']km']);
disp('~~~~~~~~');

% ~~~~ Question 7 ~~~~
disp(['Q7) The latitude of directTV10 at the initial time is [',num2str(latDms),']deg min sec'])
disp(['    The longitude of directTV10 at the initial time is [',num2str(longDms),']deg min sec'])

%% Orbital Propagation Function
% ~~~~~~~~~~~~~~~~~~~~~~~~
function dydt = orbProp(t,y,mu)
% ~~~~~~~~~~~~~~~~~~~~~~~~
JD = t/86400;
% Ttdb = (JD - 2451545.0)/36525;
R12 = y(1:3); R12n = norm(R12);
[R13] = moon(JD); R13n = norm(R13); R13 = transpose(R13);
[R14] = sun(JD)*149597900; R14n = norm(R14); R14 = transpose(R14);
R32 = R12-R13; R32n = norm(R32);
R42 = R12-R14; R42n = norm(R42);
R23 = -R32;    R23n = R32n;
R24 = -R42;    R24n = R42n;

V2 = y(4:6);
A2 = -(mu(1)/R12n^3).*R12 - (mu(3)/R32n^3).*R32 - (mu(4)/R42n^3).*R42 ;
%A2 = -(mu(1)/R12n^3).*R12 + mu(3).*(R23./R23n^3-R13./R13n^3) + mu(4).*(R24./R24n^3-R14./R14n^3);
dydt = [V2;A2];
end

%% Position of Moon and Sun 
% Personal functions for finding the position of sun and moon
% function [r13,r13n] = moon(Ttdb)
% % ~~~~ Moon Variables ~~~~
% lamEclipMoon = 218.32 + 481267.8813*Ttdb + 6.29*sind(134.9+477168.85*Ttdb)...
%                 - 1.27*sind(259.2-413335.38*Ttdb) + 0.66*sind(235.7+890534.23*Ttdb)...
%                 + 0.21*sind(269.9+954397.70*Ttdb) - 0.19*sind(357.5+35999.05*Ttdb)...
%                 - 0.11*sind(186.6+966404.05*Ttdb);
% phiEclipMoon = 5.13*sind(93.3+483202.03*Ttdb) + 0.28*sind(228.2+960400.87*Ttdb)...
%                 - 0.28*sind(318.3+6003.18*Ttdb) - 0.17*sind(217.6-407332.20*Ttdb);
% P = 0.9508 + 0.0518*cosd(134.9+477198.85*Ttdb)...
%      + 0.0095*cosd(259.2-413335.38*Ttdb) + 0.0078*cosd(235.7+890534.23*Ttdb)...
%      + 0.0028*cosd(269.9+954397.7*Ttdb);
% EMoon = 23.439291 - 0.0130042*Ttdb - 1.64E-7*Ttdb^2 + 5.04E-7*Ttdb^3;
% 
% % ~~~~ Moon Position ~~~~
% r13n = 1/sind(P); 
% r13 = r13n.*[cosd(phiEclipMoon)*cosd(lamEclipMoon);
%              cosd(EMoon)*cosd(phiEclipMoon)*sind(lamEclipMoon)-sind(EMoon)*sind(phiEclipMoon);
%              sind(EMoon)*cosd(phiEclipMoon)*sind(lamEclipMoon)+cosd(EMoon)*sind(phiEclipMoon)];
% end
% 
% function [r14,r14n] = sun(Ttdb)
% % ~~~~ Sun Variables ~~~~
% lamMSun = 280.460 + 36000.771 * Ttdb;       % (deg)
% Msun = 357.5291092 + 35999.05034 * Ttdb;    % True anomoly of sun(deg)
% lamEclipSun = lamMSun + 1.914666471 * sind(Msun) + 0.019994643 * sind(2*Msun);
% ESun = 23.439291 - 0.0130042 * Ttdb;
% 
% % ~~~~ Sun Position~~~~
% r14n = (1.000140612 - 0.016708617*cosd(Msun) - 0.000139589*cosd(2*Msun))*149597900; % (km)
% r14 = r14n.*[cosd(lamEclipSun);...
%               cosd(ESun) * sind(lamEclipSun);...
%               sind(ESun) * sind(lamEclipSun)];
% end

%% Change of Frame Functions
function c1 = C1(a)
c1 = [1, 0, 0;...
    0, cosd(a), sind(a);...
    0, -sind(a), cosd(a)];
end

function c3 = C3(a)
c3 = [cosd(a), sind(a), 0;...
    -sind(a), cosd(a), 0;...
    0, 0, 1];
end

