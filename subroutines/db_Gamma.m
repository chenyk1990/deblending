 function d=yc_Gamma(Array,rt,nt,nx,ny,ntpad,nxpad,nypad)
% yc_Gamma: Gamma operator for blending
% 
% ss = size(Array);
% ntpad = ss(1);
% nxpad=ss(2);
% nypad=ss(3);

 lengthdd = max(max(rt)) + nt;

 d = zeros(lengthdd,1);
 for iyshot=1:ny
 for ishot=1:nx
   istart = rt(ishot,iyshot) ;
   for itime=1:nt
      d(itime+istart,1) = d(itime+istart,1) + Array(itime,ishot,iyshot);
   end
 end
 end


     
