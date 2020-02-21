%parameters

I=1000;%y resolution
J=1000;%x resolution
%[I,J]=size(Z1);
delta_x=1000;%spacing
delta_y=1000;
m=.5;%m in SP law
k=1e-5;%erodibility
dt=1e5;%time step
z=rand(I,J);%initial grid
t_end=10e6;%total time
U=.0001;
%% Calculate slope
    tic;

for t=0:dt:t_end
slopes=zeros(size(z)); %slope to receiver
dxs=zeros(size(z))+999999999999;%distances to receiver
r=zeros(size(z));%receivers
for i=2:I-1
    for j=2:J-1 %loop through all nodes
        slopemax=0;
        for ip=-1:1 %loop through all neighbors
            for jp=-1:1
                slope=(z(i,j)-z(i+ip,j+jp))/sqrt((delta_y*jp).^2+(delta_x*ip).^2+1e-8);
                if slope > slopemax %Find maximum slope (>0)
                    slopes(i,j)=slope;%slopes
                    dxs(i,j)=sqrt((delta_y*jp).^2+(delta_x*ip).^2+1e-8); %distances
                    r(i,j)=(j+jp-1)*I+i+ip; %assign linear index to max slope
                    slopemax=slope;
                end
            end
        end
    end
end

%% Create stack

s=zeros(1,numel(z)); %Stack
edges=find(slopes==0); %Edge nodes
[edgei,edgej]=ind2sub(size(z),edges);
c=1;
c2=1;

for l=1:length(edgei) %Loop through all local sinks and edge nodes
    i=edgei(l);
    j=edgej(l);
    s(c2)=(j-1)*I+i;%assign edge or sink to end of stack
    c2=c2+1;

    while c<=c2&&c<numel(z) %while the current node is less than length of stack...

    for ip=-1:1
        for jp=-1:1 %Loop through all neighbors of current node
            if j+jp>1&&i+ip>1&&i+ip<I&&j+jp<J
                if r(i+ip,j+jp)==(j-1)*I+i % If neighbor is donor to current node
                    s(c2)=(j+jp-1)*I+i+ip;  %Add neighbor to stack
                    c2=c2+1;
                end
                
            end
        end
    end
    i=mod(s(c)-1,I)+1; % Go to next node in stack.
    j=floor((s(c)-1)/I)+1;
    c=c+1;
    if c>numel(z)
        break
    end

    end 
end
s(s==0)=[];


%% Calculate Drainage area
A=ones(size(z))*delta_x*delta_y;r(r==0)=1;
for l=length(s):-1:1
    A(r(s(l)))=A(r(s(l)))+A(s(l)); %cumulatively sum downstream
end
A(1)=1;
%% Calculate topography
for l=1:length(s)
    z(s(l))=1/(1+k*dt*A(s(l))^m/dxs(s(l)))*(z(s(l))...
        +k*dt*A(s(l))^m/dxs(s(l))*z(r(s(l)))); %implicit equation, goes upstream
end


z=z+U*dt;
z(end,:)=0;%Boundary conditions
z(:,1)=0;
z(:,end)=0;
z(1,:)=0;
if mod(t/dt,5)==0
figure(1);pcolor(z);shading interp;colormap jet;drawnow;
end
t
end
toc;

