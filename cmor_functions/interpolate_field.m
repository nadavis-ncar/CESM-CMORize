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
interp_len=length(dim.interp);
field_size(interp_dim)=interp_len;
field_out=zeros(field_size);

if interp_dim==2
   loop_num=3;
else
   loop_num=2;
end

%Interpolate
disp('Interpolating')
for i=1:size(field,1)
   disp([sprintf('%0.2f',100*i/size(field,1)),'% complete'])
   if do_native==1
      inner_loop=size(field,3);
      parfor j=1:size(field,loop_num)
         local_field=squeeze(field(i,j,:,:));
         local_ps=squeeze(ps(i,j,:));
         local_field_out=zeros(size(local_field,3),interp_len);
         for k=1:inner_loop
            local_field_out(k,:)=interp1(squeeze(local_ps(k))*b+a,squeeze(local_field(k,:)),...
                                         dim.interp,'linear',missing_value);
         end
         field_out(i,j,:,:)=local_field_out;
      end
   else
      fieldsize=length(size(field));
      switch fieldsize
         case 2
            parfor j=1:size(field,2)
               field_out(i,j)=interp1(dim.native.value,...
                                      squeeze(field(i,j,:)),dim.interp,...
                                      'linear',missing_value);
            end
         case 3
            parfor j=1:size(field,3)
               field_out(i,:,j)=interp1(dim.native.value,...
                                        squeeze(field(i,:,j)),dim.interp,...
                                        'linear',missing_value);
            end
         case 4
         inner_loop=size(field,3);
         parfor j=1:size(field,2)
            local_field=squeeze(field(i,j,:,:));
            local_field_out=zeros(inner_loop,interp_len);
            for k=1:inner_loop
               local_field_out(k,:)=interp1(dim.native.value,...
                                          squeeze(local_field(k,:)),dim.interp,...
                                         'linear',missing_value);
            end
            field_out(i,j,:,:)=local_field_out;
         end
      end
   end
end
