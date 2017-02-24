%% Mantle Tomography Model
%
%---------------------------------------------------------------
% This script plots the mantle tomography model                 |
% SL2013NA                                                      |
%                                                               |
% Last updated 2/16/2017 by aburky@princeton.edu                |
% Last updated 2/16/2017 by fjsimons@princeton.edu              |
%---------------------------------------------------------------

% Where is the model? User will change
filepath='/home/aburky/work/bermuda/mantleModels/SL2013NA_tri-grd_v1.0';

% Load the model into the variable, model
model=dlmread('SL2013NA_n-k.mod','',1,0);

% Create scattered interpolant of model
% x values (longitudes, deg)
x=model(:,2);
% y values (latitudes, deg)
y=model(:,3);
% z values (depths, km)
z=model(:,1);
% vs values (m/s)
vs=model(:,7);

% Create interpolant so model values can be found for any (lat,lon)
% in the domain, also convert Vs to km/s
Vs=scatteredInterpolant(x,y,z,(vs/1000));

%% Choose the starting and end points of the cross section
% (in degrees)
lat1=39;
lon1=280;
lat2=28;
lon2=302;

% Choose number of (lat,lon) points and number of depth points
nxy=1000;
nz=1000;

[lat, lon]=track2(lat1,lon1,lat2,lon2,'','degrees',nxy);
lon=lon+360;

% By default the depth ranges from 20 to 585 km
% this should not be changed since the model is only valid in
% this domain
dep=linspace(20,585,nz);

%% Construct a grid of model values
% Rows correspond to fixed depths
% Columns correspond to fixed (lat,lon)

% Pre-allocate matrix of absolute Vs
absVs=zeros(nz,nxy);

for i=1:nz
    for j=1:nxy
        absVs(i,j)=Vs(lon(j),lat(j),dep(i));
    end
end

%% Compute 1-D velocity average for the model region
 
% Pre-allocate variables
meanVs=zeros(1,nz);
delVs=zeros(nz,nxy);

% Construct array of averages for each depth slice
for i = 1:nz
    meanVs(i) = mean(absVs(i,:));
end
% Replace absolute velocities with relative velocities
for i = 1:nz
    delVs(i,:) = ((absVs(i,:)-meanVs(i))./meanVs(i))*100;
end

%% Make coordinates for plotting
arclen=distance('gc',[lat1,lon1],[lat2,lon2]);
depth=6371-dep;
x = zeros(nz,nxy);
y= zeros(nz,nxy);

% Convert polar coordinates to cartesian values
for ir = 1:nz;
    for ith = 1:nxy;
        r = depth(ir);
        th = ((pi/2)-((arclen/2)*(pi/180)))+(ith-1)*((arclen/(nxy-1))*(pi/180));
        x(ir,ith) = r*cos(th);
        y(ir,ith) = r*sin(th);
    end
end

% Rotate the velocity grid to give it the right sense (relative)
% flipVs=transpose(delVs);
rotVs=rot90(delVs,1);
flipVs=transpose(rotVs);

%% Make the tomographic cross section plot
% Load the color palette
load('kelim.mat');

% Make the cross section plot
pcolor(x,y,flipVs)
shading interp
hold on
[C,H]=contour(x,y,flipVs,7,'k');
set(H,'linewidth',0.01)
ttl=title('SL2013NA Upper Mantle Cross-Section');
P=get(ttl,'Position');
set(ttl,'Position',[P(1) P(2)+100 P(3)]);
ylabel('Depth (km)')
xlabel('Surface Distance (km)')
h = colorbar('location','southoutside');
set(h,'position',[0.27 0.12 .5 .07])
caxis([-4 4]);
ylabel(h, 'V_{S} Variation from Section 1-D Average (%)')
colormap(kelim)
axis equal tight
axis off

% Add the location of Bermuda to the plot
br=6385;
bth=((pi/2)+((arclen/2)*(pi/180)))-(14.4525*(pi/180));
bx=br*cos(bth);
by=br*sin(bth);
plot(bx,by,'k^','MarkerFaceColor','k');

% Add dots to corresponding to map dots
r1=6371;
th1=((pi/2)+((arclen/2)*(pi/180)));
x1=r1*cos(th1);
y1=r1*sin(th1);
plot(x1,y1,'ko','MarkerFaceColor','w','MarkerSize',10);

r2=6371;
th2=((pi/2)+((arclen/4)*(pi/180)));
x2=r2*cos(th2);
y2=r2*sin(th2);
plot(x2,y2,'ko','MarkerFaceColor','[0.75 0.75 0.75]','MarkerSize',10);

r3=6371;
th3=(pi/2);
x3=r3*cos(th3);
y3=r3*sin(th3);
plot(x3,y3,'ko','MarkerFaceColor','[0.5 0.5 0.5]','MarkerSize',10);

r4=6371;
th4=((pi/2)-((arclen/4)*(pi/180)));
x4=r4*cos(th4);
y4=r4*sin(th4);
plot(x4,y4,'ko','MarkerFaceColor','[0.25 0.25 0.25]','MarkerSize',10);

r5=6371;
th5=((pi/2)-((arclen/2)*(pi/180)));
x5=r5*cos(th5);
y5=r5*sin(th5);
plot(x5,y5,'ko','MarkerFaceColor','k','MarkerSize',10);

% Add a box around the figure and axes

% Pre-allocate variables
xl=zeros(2000,1);
yl=zeros(2000,1);
xr=zeros(2000,1);
yr=zeros(2000,1);
xb=zeros(1,2000);
yb=zeros(1,2000);
xt=zeros(1,2000);
yt=zeros(1,2000);
x410=zeros(1,2000);
y410=zeros(1,2000);

% Right boundary
ith=1;
for ir = 1:2000
    r = min(depth)+(ir-1)*((6371-min(depth))/1999);
    th = ((pi/2)-((arclen/2)*(pi/180)))+(ith-1)*((arclen/1999)*(pi/180));
    xr(ir,ith) = r*cos(th);
    yr(ir,ith) = r*sin(th);
end
plot(xr,yr,'k');

% Left boundary
for ir = 1:2000
    r = min(depth)+(ir-1)*((6371-min(depth))/1999);
    th = ((pi/2)+((arclen/2)*(pi/180)))+(ith-1)*((arclen/1999)*(pi/180));
    xl(ir,ith) = r*cos(th);
    yl(ir,ith) = r*sin(th);
end
plot(xl,yl,'k');

% Bottom boundary
for ir = 1;
    for ith = 1:2000
        r = min(depth);
        th = ((pi/2)-((arclen/2)*(pi/180)))+(ith-1)*((arclen/1999)*(pi/180));
        xb(ir,ith) = r*cos(th);
        yb(ir,ith) = r*sin(th);
    end
end
plot(xb,yb,'k');

% Top boundary
for ir = 1;
    for ith = 1:2000
        r = 6371;
        th = ((pi/2)-((arclen/2)*(pi/180)))+(ith-1)*((arclen/1999)*(pi/180));
        xt(ir,ith) = r*cos(th);
        yt(ir,ith) = r*sin(th);
    end
end
plot(xt,yt,'k');

% 410 Boundary
for ir = 1;
    for ith = 1:2000
        r = 5961;
        th = ((pi/2)-((arclen/2)*(pi/180)))+(ith-1)*((arclen/1999)*(pi/180));
        x410(ir,ith) = r*cos(th);
        y410(ir,ith) = r*sin(th);
    end
end
plot(x410,y410,'k--');

% Label the discontinuities
% 410 km discontinuity
r410=5970;
th410=((pi/2)+((arclen/2)*(pi/180)))+(2*(pi/180));
xtxt410=r410*cos(th410);
ytxt410=r410*sin(th410);
txt410='410';
text(xtxt410,ytxt410,txt410);

hold off;

%% Make plot of 1-D Velocity Model

figure('Position',[100, 100, 300, 600]);
plot(meanVs,dep)
set(gca,'XAxisLocation','top','ydir','reverse')
xlabel('Regional Average V_{S} (km/s)')
ylabel('Depth (km)')
axis([3 7 20 585])
