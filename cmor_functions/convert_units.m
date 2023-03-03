%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert units 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=convert_units(data,ps,a,b,conversion)

if strcmp(conversion,'number_density_to_molar_mixing_ratio')
   parfor i=1:size(data,1)
      data_local=squeeze(data(i,:,:,:));
      ps_local=squeeze(ps(i,:,:));
      data_local_out=zeros(size(data_local));
      for j=1:size(data_local,1)
         for k=1:size(data_local,2)
            for t=1:size(data_local,3)
               data_local_out(j,k,t)=(100^3)*squeeze(data_local(j,k,t))/((6.02e23/0.028)*...
                                        (a(k)+b(k)*squeeze(ps_local(j,t)))/287);
            end
         end
      end
      data(i,:,:,:)=data_local_out;
   end
elseif strcmp(conversion,'molar_mixing_ratio_to_DU')
   data=data*(287*273.15/1e5);
else
   error('unknown command passed to convert_units')
end

%Future work: use json with value/key pairs to determine conversions for a broader array of options; 
%have call option for no supplied conversion (just check). We can plan ahead 
%and ensure that var_dictionary loads temp when there is a chemical integral.
