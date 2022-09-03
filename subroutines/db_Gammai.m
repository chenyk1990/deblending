 function Array=yc_Gammai(d,rt,nt,nx,ny,ntpad,nxpad,nypad)
% yc_Gammai: adjoint Gamma operator for combing

 lengthdd = max(max(rt)) + nt;

 Array = zeros(ntpad,nxpad,nypad);
 for iyshot=1:ny
 for ishot=1:nx
   istart = rt(ishot,iyshot);
   for itime=1:nt
      Array(itime,ishot,iyshot) = Array(itime,ishot,iyshot) + d(itime+istart,1);
   end
 end
 end


     
