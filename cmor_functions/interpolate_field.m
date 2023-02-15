%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Interpolates along specified dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function field_out=interpolate_field(field,interp_dim,dim,missing_value,varargin)

do_native=0;

%Optional: load in surface pressure, a/b's for vertical interp, generate pressure array
if nargin>4
   a=varargin{1};
   b=varargin{2};
   ps=varargin{3};
   do_native=1;
end



%Set up output field
field_size=size(field);
field_size(interp_dim)=length(dim.interp);
field_out=zeros(field_size);
if do_native==1
   field_out_local=zeros(size(field_out,3),size(field_out,4));
end

if interp_dim==2
   loop_num=3;
else
   loop_num=2;
end

%Interpolate
disp('Interpolating')
for i=1:size(field,1)
   disp([sprintf('%0.2f',100*i/size(field,1)),'% complete'])
   for j=1:size(field,loop_num)
      if do_native==1
         size(field_out)
size(field)
size(field_out_local)
         if length(size(field_out))==3
            parfor k=1:size(field,3)
               field_out_local(k)=interp1(squeeze(ps(i,j,k))*b+a,field(i,j,k,:),...
                                  dim.interp,'linear',missing_value);
            end
            field_out(i,j,:)=field_out_local;
         else
            parfor k=1:size(field,3)
               field_out_local(k,:)=interp1(squeeze(ps(i,j,k))*b+a,field(i,j,k,:),...
                                  dim.interp,'linear',missing_value);
            end
         end
         field_out(i,j,:,:)=field_out_local;
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
