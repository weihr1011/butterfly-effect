clear
clf
%initialize stuff
par=MyFunc;sim=simulation(par.Ngames);dis_all=zeros(1,par.tmax/200);
nFrame = numel(1:200:par.tmax);frame(nFrame)= struct('cdata',[],'colormap',[]);
fig=inital_fig;
%start doing calculation and plot
for t=1:par.tmax%loop over all time steps
    cla(fig.ax2)
    dis=zeros(2,100,25);%initialize a 2x100x25 array for dispersion
    for n=1:par.Ngames%loop over all games
        hold on
        sim(n)=verlet(sim(n));%update position verlet,velocity verlet       
        dis(:,:,n)=sim(n).x;%store all x for n games to calculate dispersion
        if mod(t,200)==0%plot every 200 steps
            subplot(fig.ax2)            
            plot(fig.ax2,sim(n).x(1,:),sim(n).x(2,:),'.','MarkerSize',10);
            frame(t)=getframe(gcf);
            if n==25%calculate dispersion after getting all 25 games x info
                disper=dispersion(dis);
                dis_all(t/200)=disper;
                subplot(fig.ax1)                              
                addpoints(fig.h,t*1e-5,dis_all(t/200))
                drawnow
            end           
        end
    end                               
end
dis_all=log(dis_all);
%Polyfit(dis_all)
getvideo(frame)%generate mp4

function par=boundaryChecking(par)%check if exceed boundary
%check if exceed x=1
dx1=par.x(1,:)-par.xmax(1,:);
indexbounce=dx1>0;
par.x(1,indexbounce)=par.xmax(1,indexbounce)-dx1(indexbounce);
par.v(1,indexbounce)=-par.v(1,indexbounce);
%check if exceed x=0
dx2=par.x(1,:)-par.xmin(1,:);
indexbounce=dx2<0;
par.x(1,indexbounce)=par.xmin(1,indexbounce)-dx2(indexbounce);
par.v(1,indexbounce)=-par.v(1,indexbounce);
%check if exceed y=1
dy1=par.x(2,:)-par.xmax(2,:);
indexbounce=dy1>0;
par.x(2,indexbounce)=par.xmax(2,indexbounce)-dy1(indexbounce);
par.v(2,indexbounce)=-par.v(2,indexbounce);
%check if exceed y=0
dy2=par.x(2,:)-par.xmin(2,:);
indexbounce=dy2<0;
par.x(2,indexbounce)=par.xmin(2,indexbounce)-dy2(indexbounce);
par.v(2,indexbounce)=-par.v(2,indexbounce);
end

function par=MyFunc %initiate variables
fixrng
par.N=100; %numbers of particles
par.tmax=1e5; %max time
par.dt=1e-5; %time step
par.x=rand(2,par.N); %inital position
par.v=randn(2,par.N);%initial velocity
par.xmax=ones(2,par.N)-1e-2; %max x,y
par.xmin=zeros(2,par.N)+1e-2;%min x,y
par.Ngames=25; %number of games
par.size=10;%plot size
par.eps=1e-12;%gaussian random pertubation
par.k=zeros(2,100); % Function handle for forces.
par.alpha=2; % Force power law.
par.F0=20; % Force kernel magnitude.
par.r0=0.1; % Force core radius (like a particle size)
par.Fx=zeros(1,100);%initialize Fx
par.Fy=zeros(1,100);%initialize Fy
par.prev_Fx=zeros(1,100);%initialize Fx from last step
par.prev_Fy=zeros(1,100);%initialize Fy from last step
end
function sim=simulation(n)%iteration of n games
par=MyFunc;
for i=1:n
    if i==1
        sim(i)=par;
    else
        sim(i)=par;
        sim(i).v=par.v+rand(1)*par.eps;
    end
end
end
function fixrng %generating same random value each time
%https://www.mathworks.com/help/matlab/ref/rng.html
rng('default')
end
function d=distance(par)%calculate distance between each particle
d.x=par.x(1,:)'-par.x(1,:);
d.y=par.x(2,:)'-par.x(2,:);
d.d=sqrt(d.x.^2+d.y.^2);
end
function F=inv_rsq(sim)% softened finite range force law
d=distance(sim);
Fx=sim.F0.*max(0,1./(sim.r0.^sim.alpha+d.d.^sim.alpha)-1/(2*sim.r0^sim.alpha)).*d.x./(sim.r0+d.d);
Fy=sim.F0.*max(0,1./(sim.r0.^sim.alpha+d.d.^sim.alpha)-1/(2*sim.r0^sim.alpha)).*d.y./(sim.r0+d.d);
Fx=sum(Fx,2);
Fy=sum(Fy,2);
F.x=Fx';
F.y=Fy';
end
function F=inv_rsq_2(sim)% softened 'infinite range'force law
d=distance(sim);
Fx=sim.F0./(sim.r0.^sim.alpha+d.d.^sim.alpha).*d.x./(sim.r0.^sim.alpha+d.d);
Fy=sim.F0./(sim.r0.^sim.alpha+d.d.^sim.alpha).*d.y./(sim.r0.^sim.alpha+d.d);
Fx=sum(Fx,2);
Fy=sum(Fy,2);
F.x=Fx';
F.y=Fy';
end
function F=inv_rsq_3(sim)% softened finite range force law
d=distance(sim);
Fx=sim.F0.*max(0,1./(d.d.^sim.alpha)-1/(sim.r0^sim.alpha)).*d.x./(sim.r0+d.d);
Fy=sim.F0.*max(0,1./(d.d.^sim.alpha)-1/(sim.r0^sim.alpha)).*d.y./(sim.r0+d.d);
Fx=sum(Fx,2);
Fy=sum(Fy,2);
F.x=Fx';
F.y=Fy';
end
function par=verlet(par)%velocity&position vetlet
par.x(1,:)=par.x(1,:)+par.v(1,:)*par.dt+1/2.*(par.Fx.*par.dt^2);%update X
par.x(2,:)=par.x(2,:)+par.v(2,:)*par.dt+1/2.*(par.Fy.*par.dt^2);%update Y
par=boundaryChecking(par);
K=inv_rsq(par);%get force
par.Fx=K.x;%Fx
par.Fy=K.y;%Fy
par.v(1,:)=par.v(1,:)+1/2*(par.Fx+par.prev_Fx)*par.dt;%update Vx
par.v(2,:)=par.v(2,:)+1/2*(par.Fy+par.prev_Fy)*par.dt;%update Vy
par.prev_Fx=K.x;%stored Fx for next time step calculation
par.prev_Fy=K.y;%stored Fy for next time step calculation
end
function dis=dispersion(par)%calculate dispersion
ave=sum(par,3)./25;%calculate average of x and y for each ball(2x100 array)
ave=par-ave;%calculate x-x_bar,y-y_bar(2x100x25 array)
ave=ave.^2;%(x-x_bar)^2,(y-y_bar)^2(2x100x25 array)
dis=sum(ave,3)./25;%2x100 array, first row is variance of x, second row is vriance of y
dis=sum(dis,2);%first row is total variance of x,second row is total variance of y
dis=sum(dis);% dispersion square
dis=sqrt(dis);%dispersion
end
function getvideo(par)%generate mp4 video
myWriter=VideoWriter('butterfly4','MPEG-4');
myWriter.FrameRate=20;
open(myWriter);
writeVideo(myWriter,par(200:200:end));
close(myWriter)
end
function fig=inital_fig%initialize plot
fig.ax1=subplot(1,2,1);
xlabel('time')
ylabel('dispersion')
set(gca, 'YScale', 'log')
fig.h=animatedline('Marker','.');
axis square
title('dispersion')
fig.ax2=subplot(1,2,2);
axis square
axis([0 1 0 1])
title('particles')
end
function Polyfit(par)
hold off
subplot(1,2,1)
t=linspace(0,1,5);
p=polyfit(t,par,1);
y1=polyval(p,t);
hold on
plot(t,par,'.')
plot(t,y1)
hold off
end
