%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Expand surface pressure to full pressure field with hybrid coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pres=calculate_pressure(ps,a,b,varargin)

if nargin>3
   column=varargin{1};
else
   column=0;
end

if column==0
   pres=zeros(size(ps,1),size(ps,2),length(a),size(ps,3));
   parfor i=1:size(pres,1)
      ps_local=squeeze(ps(i,:,:));
      pres_local=zeros(size(ps_local,1),size(ps_local,2),length(a));
      for j=1:size(ps_local,1)
         for t=1:size(ps_local,2)
            pres_local(j,t,:)=a(:)+b(:)*squeeze(ps_local(j,t));
         end
      end
      pres(i,:,:,:)=permute(pres_local,[1 3 2]);
   end
else
   pres=a(:)+b(:)*ps;
end
