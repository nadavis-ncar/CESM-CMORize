%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolates along specified dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field_out=interpolate_field_out(field,interp_dim,dim,missing_value,varargin)

disp('in func')
do_native=0;
%Optional: load in surface pressure, a/b's for vertical interp, generate pressure array
if nargin>4
   a=varargin{1};
   b=varargin{2};
   ps=varargin{3};
   do_native=1;
end
disp('setting up')
%Set up output field
field_size=size(field);
field_size(interp_dim)=length(dim.interp);
field_out=zeros(field_size);

if do_native==1
   local_field_out=zeros(size(field_out,4),size(field_out,3));
end

if interp_dim==2
   loop_num=3;
else
   loop_num=2;
end

disp('starting')
%Interpolate
for i=1:size(field,1)
   disp([sprintf('%0.2f',100*i/size(field,1)),'% complete'])
   for j=1:size(field,loop_num)
      if do_native==1
         local_field=squeeze(field(i,j,:,:))';
         local_dim=squeeze(dim.native.value(i,j,:,:))';
         parfor k=1:size(field,4)
           local_field_out(k,:)=interp1(squeeze(ps(i,j,k))*b+a,...
                                        squeeze(field(i,j,:,k),...
                                        dim.interp,'linear',missing_value);
         end
         field_out(i,j,:,:)=local_field_out';
     else
         fieldsize=length(size(field_out));
         switch fieldsize
            case 2
               field_out(i,j)=interp1(dim.native.value,...
                                  squeeze(field(i,j,:)),dim.interp,...
                                  'linear',missing_value);
            case 3
               field_out(i,:,j)=interp1(dim.native.value,...
                                  squeeze(field(i,:,j)),dim.interp,...
                                  'linear',missing_value);
            case 4
               for k=1:size(field_out,3)
                  field_out(i,j,k,:)=interp1(dim.native.value,...
                                  squeeze(field(i,j,k,:)),dim.interp,...
                                  'linear',missing_value);
               end
         end
      end
   end
end

