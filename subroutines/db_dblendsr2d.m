function [dout,parm2] = db_dblendsr2d(din,mode,parm);
% 2D deblending using fk thresholding ( in the source and reciever domain)
%  IN   	in:  intput combed blended data
%
% 		parm: input parameters
%  		mode=1,	din: long node trace (ntmax,nr), no unblended data
%			2,	din: long node trace (ntmax,nr), with unblended data (combed) (nt,nx,nr), output snr(niter,nr)
%			3,	din: combed data (nt,nx,nr), no unblended data
%			4,	din: combed data (nt,nx,nr), with unblended data (combed) (nt,nx,nr), output snr(niter,nr)
%			5,   din: unblended data (nt,nx,nr), output is long node trace
%			(ntmax,nx,nr), without maxfold
%			6,   din: unblended data (nt,nx,nr), output is combed data
%			(nt,nx,nr), with maxfold
%			7,   din: unblended data (nt,nx,nr), output is combed data (nt,nx,nr),
%           
%  OUT   out:  output deblended data
%
%  Copyright (C) 2015 The University of Texas at Austin
%  Copyright (C) 2015 Yangkang Chen
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
if mode <=4
    nt=parm.nt;
    nx=parm.nx;
    nr=parm.nr;
    rt=parm.rt;
    niter=parm.niter;
    a=parm.a;
    alpha=parm.alpha;
    lambda=parm.lambda;
    ntpad=parm.ntpad;
    nxpad=parm.nxpad;
    verb=parm.verb;
    %    lambda=0.7 should be given
    if mode==2 || mode==4
        ub=parm.ub;
        snr=zeros(niter,nr);
    else
        mis=zeros(niter,nr);
    end
end

if mode >=5
    nt=parm.nt;
    nx=parm.nx;
    nr=parm.nr;
    ntpad=parm.ntpad;
    nxpad=parm.nxpad;
    if mode~=7 
        maxfold=parm.maxfold;
        seed=parm.seed;
    end
    verb=parm.verb; 
    if mode==7 
        rt=parm.rt;
    end
end

% default setting
ny=1;
nypad=ny;

% allocate memory
m_n = zeros(nt,nx) ;
d_b = zeros(nt,nx) ;
d_res = zeros(nt,nx) ;
mprime = zeros(nt,nx) ;
m_up = zeros(nt,nx) ;

if mode<=4
    % calculate the maximum blending fold
    B = ones(ntpad,nxpad,nypad) ;
    data = db_Gamma(B,rt,nt,nx,ny,ntpad,nxpad,nypad);
    maxfold = max(max(data));
    
    if verb==1
        fprintf('Maxfold is %d \n', maxfold);
    end
end

if mode<=4
    
    
    
    dout=zeros(nt,nx,nr);
    for ir=1:nr
        
        if mode>=3
            d_b_com=din(:,:,ir);
            %         vmx=prctile(reshape(abs(fftn(d_b_com)),nt*nx,1),99.999999);
            vmx=max(max(abs(fftn(d_b_com))));
            m_n1=zeros(size(d_b_com));
            m_up=zeros(size(d_b_com));
        else
            d_b_com=db_Gammai(din(:,ir),rt,nt,nx,1,ntpad,nxpad,1);
            %         vmx=prctile(reshape(abs(fftn(d_b_com)),nt*nx,1),99.999999);
            vmx=max(max(abs(fftn(d_b_com))));
            m_n1=zeros(size(d_b_com));
            m_up=zeros(size(d_b_com));
        end
        
        for iter=1:niter
            %%% compute the data misfit
            %d_res = d_b - d_mod;
            
            d_mod = db_Gamma(m_n1,rt,nt,nx,1,ntpad,nxpad,1);
            
            if mode<=2
                m_up=db_Gammai(din(:,ir)-d_mod,rt,nt,nx,1,ntpad,nxpad,1);
            else
                d_mod_com = db_Gammai(d_mod,rt,nt,nx,1,ntpad,nxpad,1);
                m_up=d_b_com-d_mod_com;
            end
            
            m_n = m_n1; % m_n and m_{n+1}
            m_n1 = m_n + m_up * lambda;
            mprime  = fftn(m_n1);
            vlim = ((niter-iter)/(niter-1)) * vmx ;
            mprime=(abs(mprime)>vlim).*mprime;
            m_n = ifftn(mprime);
            %%% m_{n+1} = (1-alpha)*m_n + alpha*A^{-1}TA(m_n+\lamda*F^{pseudoinverse}(d_obs-Fm_n))
            m_temp = m_n1 - m_n;
            m_n1 = m_n + m_temp * alpha;  % alpha=0.75/maxfold
            if mode==2 || mode==4
                %%% compute SNRs
                snr(iter,ir)=db_snr(ub(1:nt,:,ir),m_n1(1:nt,:));
                if verb==1
                    fprintf('iter = %d/%d, receiver=%d/%d\n',iter,niter,ir,nr);
                    fprintf('Current SNR = %g \n\n',snr(iter,ir));
                end
            else
                mis(iter,ir)=db_snr(m_n1(1:nt,:),m_n1(1:nt,:)+m_n(1:nt,:));
                if verb==1
                    fprintf('iter = %d/%d, receiver=%d/%d\n',iter,niter,ir,nr);
                    fprintf('Current model misfit = %g \n',mis(iter,ir));
                end
                
            end
        end
        dout(:,:,ir)=m_n1;
    end
    
    if mode==2 || mode==4
        parm2.snr=snr;
    else
        parm2.mis=mis;
    end
    
else
    if(mode==6 || mode==5)
        dfr=0.1;
        rt=gstimes(maxfold,nt,nx,ny,ntpad,nxpad,nypad,seed,dfr,1);
        % save the shot schedule
        parm2.rt=rt;
    end
    
    if mode==5
        dout=[];
        for ir=1:nr
            dtmp= db_Gamma(din,rt,nt,nx,ny,ntpad,nxpad,nypad);
            dout=[dout, dtmp];
        end
    else
        dout=zeros(ntpad,nxpad,nr);
        for ir=1:nr
            dtmp= db_Gamma(din(:,:,ir),rt,nt,nx,ny,ntpad,nxpad,nypad);
            dout(:,:,ir)=db_Gammai(dtmp,rt,nt,nx,ny,ntpad,nxpad,nypad);
        end
    end
end


