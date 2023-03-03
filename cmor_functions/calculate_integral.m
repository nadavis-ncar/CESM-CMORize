%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate integral in pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_integral(data,ps,a,b,varargin)

data_out=zeros(size(data,1),size(data,2),size(data,4));
do_trop=0;

%Tropospheric integral?
if nargin>5
   do_trop=1;
   trop=varargin{1};
end

%Loop over all columns and integrate
inner_loop_1=size(data,2);
inner_loop_2=size(data,4);

if do_trop==1
   parfor i=1:size(ps,1)
      ps_local=squeeze(ps(i,:,:));
      data_local=squeeze(data(i,:,:,:));
      data_out_local=zeros(size(data_local,1),size(data_local,3));
      trop_local=squeeze(trop(i,:,:));
      for j=1:size(data_local,1)
         for t=1:size(data_local,3)
            pres=a(:)+b(:)*squeeze(ps_local(j,t));
            [val,trop_ind]=min(abs(squeeze(trop_local(j,t)) - pres));
            if pres(trop_ind)>squeeze(trop_local(j,t))
               trop_ind=trop_ind-1;
            end
            data_out_local(j,t)=trapz(pres(trop_ind+1:end),squeeze(data_local(j,trop_ind+1:end,t)));
            data_out_local(j,t)=squeeze(data_out_local(j,t))+(pres(trop_ind+1)-squeeze(trop_local(j,t)))*squeeze(data_local(j,trop_ind,t));
         end
      end
      data_out(i,:,:)=data_out_local;
   end
else
   parfor i=1:size(ps,1)
      ps_local=squeeze(ps(i,:,:));
      data_local=squeeze(data(i,:,:,:));
      data_out_local=zeros(size(data_local,1),size(data_local,3));
      for j=1:size(data_local,1)
         for t=1:size(data_local,3)
            pres=a(:)+b(:)*squeeze(ps_local(j,t));
            data_out_local(j,t)=trapz(pres,squeeze(data_local(j,:,t)));
         end
      end
      data_out(i,:,:)=data_out_local;
   end
end
data=data_out/9.81;
