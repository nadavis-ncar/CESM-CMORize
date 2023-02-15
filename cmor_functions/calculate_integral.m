%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate integral in pressure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=calculate_integral(data,varargin)

data_out=zeros(size(data.native.value,1),size(data.native.value,2),size(data.native.value,4));

do_trop=0;

ps=varargin{1};
a=varargin{2};
b=varargin{3};

pres=calculate_pressure(ps,a,b);

%Tropospheric integral?
if nargin>5
   if strcmp('tropopause',varargin{4})
      do_trop=1;
      trop=data.native_secondary.value;
   end
end

%Loop over all columns and integrate
for i=1:size(ps,1)
   disp([sprintf('%3.0f',100*i/size(ps,1)),'% complete'])
   for j=1:size(ps,2)
      parfor t=1:size(data.native.value,4)
         if do_trop==1
            [val,trop_ind]=min(abs(squeeze(trop(i,j,t)) - squeeze(pres(i,j,t,:))));
            if pres(i,j,t,trop_ind)>trop(i,j,t)
               trop_ind=trop_ind-1;
            end
            data_out(i,j,t)=trapz(squeeze(pres(i,j,trop_ind+1:end,t)),data.native.value(i,j,trop_ind+1:end,t));
            data_out(i,j,t)=data_out(i,j,t)+(pres(i,j,trop_ind+1,t)-trop(i,j,t))*data.native.value(i,j,trop_ind,t);
         else
            data_out(i,j,t)=trapz(squeeze(pres(i,j,:,t)),squeeze(data.native.value(i,j,:,t)));
         end
      end
   end
end

data.native.value=data_out/9.81;
