%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert omega to w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data=omega_to_w(data,ps,a,b)

for i=1:size(data,1)
   for j=1:size(data,2)
      for t=1:size(data,4)
         data(i,j,:,t)=(9.81*(a(:)+b(:)*squeeze(ps(i,j,t)))/287).*squeeze(data(i,j,:,t));
      end
   end
end
