function [ rt ] = gstimes(maxfold,nt,nx,ny,ntpad,nxpad,nypad,seed,dfr,mode)
% generate shottime for prespecified max blending fold
% when maxfold >=3, the ISS is thought to has a very high blending fold.
% mode=1:	ascending
%     =2:   random

if(nargin<7)
    error('Please input at least 6 parameters!');
end
if(nargin==7)
    seed=201314;
    dfr=.1;
    mode=2;
end

if(nargin==8)
    dfr=.1
	mode=1;
end

if(nargin==9)
	mode=1;
end

rand('state',seed);
r0=rand(nx,ny);r0=r0/max(max(r0));

if mode==1
	r0=reshape(sort(reshape(r0,nx*ny,1),'ascend'),nx,ny);
end

fr=0;
maxblendfold=-inf;
while(maxblendfold~=maxfold && maxfold>=3)
    fr=fr+dfr;
    r = r0 * nt*nx*ny *fr;
    rt = floor(r);
    
    % find max blending fold
    B = ones(ntpad,nxpad,nypad) ;
    data = yc_Gamma(B,rt,nt,nx,ny,ntpad,nxpad,nypad);
    maxblendfold = max(max(data));
    fprintf('The controling factor fr=%g (remember this value)!\n',fr);
end

if(maxfold==2)
    fr=5;
    r = r0 * nt*nx*ny *fr;
    rt = floor(r);
    B = ones(ntpad,nxpad,nypad) ;
    data = yc_Gamma(B,rt,nt,nx,ny,ntpad,nxpad,nypad);
    maxblendfold = max(max(data));
    fprintf('Max blending fold is %d \n',maxblendfold);   
end

if(maxfold==1)
    fr=10;
    r = r0 * nt*nx*ny *fr;
    rt = floor(r);
    B = ones(ntpad,nxpad,nypad) ;
    data = Gamma1(B,rt,nt,nx,ny,ntpad,nxpad,nypad);
    maxblendfold = max(max(data));
    fprintf('Max blending fold is %d \n',maxblendfold);   
end
end

