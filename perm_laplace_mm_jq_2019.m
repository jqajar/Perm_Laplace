% Permeability estimation from 3D µ-CT core images
% based on the solution of the Laplace's equation for pressure

% Written by: Ms. Maryam Mohammadi and Dr. Jafar Qajar

% Copyright 2019

% This function receives a binarized image as "inputmedia" and
% calculates permeability of the medium.

function KK = perm_laplace(inputmedia)

[x y z]=size(inputmedia);
Pstart=200;
Pfinish=100;
containerr=bwconncomp(1-inputmedia,6);
container=labelmatrix(containerr);
nodemedia=zeros(x,y,z); %#ok<NASGU>
nodeno=1; % due to appropriate node number to each connected pore voxel
%%

connectedmedia=ones(x*y*z,1);
container1=zeros(x,y,z); container2=zeros(x,y,z); 
container1(:,:,1)=container(:,:,1);
container2(:,:,z)=container(:,:,z);
indx1=find(container1~=0);
indx2=find(container2~=0);
bbd=1;
bb(bbd)=10000000;
for i=1:size(indx1)
    for j=1:size(indx2)
        if container(indx1(i))==container(indx2(j))
            if  size(find(bb==container(indx1(i))),2)==0
                bb(bbd)=container(indx1(i));
                w=find(container==container(indx1(i)));% finding the pixle that connect the ** top of the media to down ** to flowing the fluid
                connectedmedia(w,1)=0;  % Number of that pixle
                bbd=bbd+1;
            end
        else      
        end
    end
end
connectmedia=reshape(connectedmedia,[x,y,z]);
%%
index=find(connectmedia==0);
[i ,j ,k]=ind2sub(size(connectmedia),index);
for ii=1:length(i)
    nodemedia(i(ii),j(ii),k(ii))=nodeno;
    nodeno=nodeno+1;
end



totalporosity=length(find(inputmedia==0))/(x*y*z);
connectporosity=length(index)/(x*y*z);
coefequ=sparse(length(index),length(index));
fixequ=sparse(length(index),1);
%%

temp=ones(size(connectmedia));
temp(2:end-1,2:end-1,2:end-1)=connectmedia(2:end-1,2:end-1,2:end-1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:));
end


temp=ones(size(connectmedia));
temp(1,2:end-1,2:end-1)=connectmedia(1,2:end-1,2:end-1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:));
end

temp=ones(size(connectmedia));
temp(x,2:end-1,2:end-1)=connectmedia(x,2:end-1,2:end-1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:));
end

temp=ones(size(connectmedia));
temp(2:end-1,1,2:end-1)=connectmedia(2:end-1,1,2:end-1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:));
end

temp=ones(size(connectmedia));
temp(2:end-1,y,2:end-1)=connectmedia(2:end-1,y,2:end-1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:));
end

temp=ones(size(connectmedia));
temp(2:end-1,2:end-1,1)=connectmedia(2:end-1,2:end-1,1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pstart;
end

temp=ones(size(connectmedia));
temp(2:end-1,2:end-1,z)=connectmedia(2:end-1,2:end-1,z);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pfinish;
end

temp=ones(size(connectmedia));
temp(1,2:end-1,1)=connectmedia(1,2:end-1,1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pstart;
end

temp=ones(size(connectmedia));
temp(x,2:end-1,1)=connectmedia(x,2:end-1,1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pstart;
end

temp=ones(size(connectmedia));
temp(2:end-1,1,1)=connectmedia(2:end-1,1,1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pstart;
end

temp=ones(size(connectmedia));
temp(2:end-1,y,1)=connectmedia(2:end-1,y,1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pstart;
end

temp=ones(size(connectmedia));
temp(1,1,2:end-1)=connectmedia(1,1,2:end-1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:));
end


temp=ones(size(connectmedia));
temp(1,y,2:end-1)=connectmedia(1,y,2:end-1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:));
end


temp=ones(size(connectmedia));
temp(x,1,2:end-1)=connectmedia(x,1,2:end-1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:));
end


temp=ones(size(connectmedia));
temp(x,y,2:end-1)=connectmedia(x,y,2:end-1);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    if connectmedia(i(m),j(m),k(m)+1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)+1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:));
end

temp=ones(size(connectmedia));
temp(1,2:end-1,z)=connectmedia(1,2:end-1,z);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pfinish;
end


temp=ones(size(connectmedia));
temp(x,2:end-1,z)=connectmedia(x,2:end-1,z);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pfinish;
end

temp=ones(size(connectmedia));
temp(2:end-1,1,z)=connectmedia(2:end-1,1,z);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)+1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)+1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pfinish;
end


temp=ones(size(connectmedia));
temp(2:end-1,y,z)=connectmedia(2:end-1,y,z);
indx=find(temp==0);
[i,j,k]=ind2sub(size(connectmedia),indx);
for m=1:length(i)
    if connectmedia(i(m)-1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)-1,j(m),k(m)))=1;
    end
    if connectmedia(i(m)+1,j(m),k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m)+1,j(m),k(m)))=1;
    end
    if connectmedia(i(m),j(m)-1,k(m))==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m)-1,k(m)))=1;
    end
    if connectmedia(i(m),j(m),k(m)-1)==0
        coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)-1))=1;
    end
    coefequ(nodemedia(i(m),j(m),k(m)),nodemedia(i(m),j(m),k(m)))=-sum(coefequ(nodemedia(i(m),j(m),k(m)),:))-1;
    fixequ(nodemedia(i(m),j(m),k(m)))=Pfinish;
end


if connectmedia(1,1,1)==0
    
    if connectmedia(2,1,1)==0
        coefequ(nodemedia(1,1,1),nodemedia(2,1,1))=1;
    end
    if connectmedia(1,2,1)==0
        coefequ(nodemedia(1,1,1),nodemedia(1,2,1))=1;
    end
    if connectmedia(1,1,2)==0
        coefequ(nodemedia(1,1,1),nodemedia(1,1,2))=1;
    end
    coefequ(nodemedia(1,1,1),nodemedia(1,1,1))=-sum(coefequ(nodemedia(1,1,1),:))-1;
    fixequ(nodemedia(1,1,1))=Pstart;
end

if connectmedia(1,1,z)==0
    if connectmedia(2,1,z)==0
        coefequ(nodemedia(1,1,z),nodemedia(2,1,z))=1;
    end
    if connectmedia(1,2,k)==0
        coefequ(nodemedia(1,1,z),nodemedia(1,2,z))=1;
    end
    if connectmedia(1,1,z-1)==0
        coefequ(nodemedia(1,1,z),nodemedia(1,1,z-1))=1;
    end
    coefequ(nodemedia(1,1,z),nodemedia(1,1,z))=-sum(coefequ(nodemedia(1,1,z),:))-1;
    fixequ(nodemedia(1,1,z))=Pfinish;
end

if connectmedia(1,y,1)==0
    if connectmedia(2,y,1)==0
        coefequ(nodemedia(1,y,1),nodemedia(2,y,1))=1;
    end
    if connectmedia(1,y-1,1)==0
        coefequ(nodemedia(1,y,1),nodemedia(1,y-1,1))=1;
    end
    if connectmedia(1,y,2)==0
        coefequ(nodemedia(1,y,1),nodemedia(1,y,2))=1;
    end
    coefequ(nodemedia(1,y,1),nodemedia(1,y,1))=-sum(coefequ(nodemedia(1,y,1),:))-1;
    fixequ(nodemedia(1,y,1))=Pstart;
end

if connectmedia(1,y,z)==0
    if connectmedia(2,y,z)==0
        coefequ(nodemedia(1,y,z),nodemedia(2,y,z))=1;
    end
    if connectmedia(1,y-1,z)==0
        coefequ(nodemedia(1,y,z),nodemedia(1,y-1,z))=1;
    end
    if connectmedia(1,y,z-1)==0
        coefequ(nodemedia(1,y,z),nodemedia(1,y,z-1))=1;
    end
    coefequ(nodemedia(1,y,z),nodemedia(1,y,z))=-sum(coefequ(nodemedia(1,y,z),:))-1;
    fixequ(nodemedia(1,y,z))=Pfinish;
end

if connectmedia(x,1,1)==0
    if connectmedia(x-1,1,1)==0
        coefequ(nodemedia(x,1,1),nodemedia(x-1,1,1))=1;
    end
    if connectmedia(x,2,1)==0
        coefequ(nodemedia(x,1,1),nodemedia(x,2,1))=1;
    end
    if connectmedia(x,1,2)==0
        coefequ(nodemedia(x,1,1),nodemedia(x,1,2))=1;
    end
    coefequ(nodemedia(x,1,1),nodemedia(x,1,1))=-sum(coefequ(nodemedia(x,1,1),:))-1;
    fixequ(nodemedia(x,1,1))=Pstart;
end

if connectmedia(x,1,z)==0
    if connectmedia(x-1,1,z)==0
        coefequ(nodemedia(x,1,z),nodemedia(x-1,1,z))=1;
    end
    if connectmedia(x,2,z)==0
        coefequ(nodemedia(x,1,z),nodemedia(x,2,z))=1;
    end
    if connectmedia(x,1,z-1)==0
        coefequ(nodemedia(x,1,z),nodemedia(x,1,z-1))=1;
    end
    coefequ(nodemedia(x,1,z),nodemedia(x,1,z))=-sum(coefequ(nodemedia(x,1,z),:))-1;
    fixequ(nodemedia(x,1,z))=Pfinish;
end

if connectmedia(x,y,1)==0
    if connectmedia(x-1,y,1)==0
        coefequ(nodemedia(x,y,1),nodemedia(x-1,y,1))=1;
    end
    if connectmedia(x,y-1,1)==0
        coefequ(nodemedia(x,y,1),nodemedia(x,y-1,1))=1;
    end
    if connectmedia(x,y,2)==0
        coefequ(nodemedia(x,y,1),nodemedia(x,y,2))=1;
    end
    coefequ(nodemedia(x,y,1),nodemedia(x,y,1))=-sum(coefequ(nodemedia(x,y,1),:))-1;
    fixequ(nodemedia(x,y,1))=Pstart;
end

if connectmedia(x,y,z)==0
    if connectmedia(x-1,y,z)==0
        coefequ(nodemedia(x,y,z),nodemedia(x-1,y,z))=1;
    end
    if connectmedia(x,y-1,z)==0
        coefequ(nodemedia(x,y,z),nodemedia(x,y-1,z))=1;
    end
    if connectmedia(x,y,z-1)==0
        coefequ(nodemedia(x,y,z),nodemedia(x,y,z-1))=1;
    end
    coefequ(nodemedia(x,y,z),nodemedia(x,y,z))=-sum(coefequ(nodemedia(x,y,z),:))-1;
    fixequ(nodemedia(x,y,z))=Pfinish;
end
presslist=-coefequ\fixequ;

pressmedia=zeros(x,y,z);
indx=find(nodemedia~=0);
[i,j,k]=ind2sub(size(nodemedia),indx);

for ii=1:length(i)
    pressmedia(i(ii),j(ii),k(ii))=presslist(ii);
end


P=pressmedia;

% pixel=2.64*1e-6;  %length of each pixel in meter
pixel=1*1e-6; %[=]m
EDT_map=bwdist(connectmedia);
k_EDT = zeros(x,y,z);
indx=find(EDT_map~=0);
for i=1:length(indx)
    k_EDT(indx(i))=((EDT_map(indx(i)) - (1/2))*pixel).^2;  % Arns and Adler 2012
end

            
A_voxel=pixel*pixel;       % Area of each pixel in squrt meter
miu=0.1;                   % viscosity of fluid in pa.s;
delX=pixel;
L=pixel*z;                 % length of media in flow direction(up to down)

%%

for i=1:x
    for j=1:y
        for k=2:z-1
            if P(i,j,k)~=0 && P(i,j,k+1)~=0 && P(i,j,k-1)~=0
                q_EDT(i,j,k)=(k_EDT(i,j,k)*A_voxel*(P(i,j,k-1)-P(i,j,k+1)))/(miu*delX*2);
                delp(i,j,k)=abs(P(i,j,k-1)-P(i,j,k+1))/2;
            end
        end
    end
end

for i=1:x
    for j=1:y
        if P(i,j,1)~=0 && P(i,j,2)~=0
            q_EDT(i,j,1)=(k_EDT(i,j,1)*A_voxel*(Pstart-P(i,j,2)))/(miu*delX*2);
            delp(i,j,1)=abs(P(i,j,1)-P(i,j,2));
        end
    end
end

for i=1:x
    for j=1:y
        if P(i,j,z)~=0 && P(i,j,z-1)~=0
            q_EDT(i,j,z)=(k_EDT(i,j,z)*A_voxel*(P(i,j,z-1)-Pfinish))/(miu*delX*2);
            delp(i,j,z)=abs(P(i,j,z-1)-P(i,j,z));
        end
    end
end

for k=1:z
    Q(k)=sum(sum(abs(q_EDT(:,:,k))));
    delP(k)=sum(sum(abs(delp(:,:,k))))/(length(find(delp(:,:,k)~=0)));
    kedtL(k)=Q(k)*miu*pixel/(A_voxel*x*y*delP(k))*1.01325e12;
    KEDTinv(k)=Q(k)*pixel/kedtL(k);
    QL(k)=Q(k)*pixel;
end
K_EDT_ser=sum(QL)/(abs(sum(KEDTinv)));

KK=sum(Q)*miu*L/(A_voxel*x*y*(Pstart-Pfinish))*1.01325e12/z;


