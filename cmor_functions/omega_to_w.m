%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert omega to w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=omega_to_w(data,ps,a,b)

parfor i=1:size(data,1)
   data_local=squeeze(data(i,:,:,:));
   ps_local=squeeze(ps(i,:,:));
   data_out_local=zeros(size(data_local));
   for j=1:size(data_local,1)
      for t=1:size(data_local,3)
         data_out_local(j,:,t)=(9.81*(a(:)'+b(:)'*squeeze(ps_local(j,t)))/287).*squeeze(data_local(j,:,t));
      end
   end
   data(i,:,:,:)=data_out_local;
end
