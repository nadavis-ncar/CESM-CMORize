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

%Interpolate
for i=1:size(field,1)

   %Data is in hybrid-sigma coordinates
   if do_native==1
      inner_loop=size(field,3);
      parfor j=1:size(field,2)
         local_field=squeeze(field(i,j,:,:));
         local_ps=squeeze(ps(i,j,:));
         if interp_dim==4
            local_field_out=zeros(inner_loop,interp_len);
            for k=1:inner_loop
               local_field_out(k,:)=interp1(local_ps(k)*b+a,squeeze(local_field(k,:)),...
                                            dim.interp,'linear',missing_value);
            end
         else
            local_field_out=zeros(interp_len,inner_loop);
            for k=1:inner_loop
               local_field_out(:,k)=interp1(local_ps(k)*b+a,squeeze(local_field(:,k)),...
                                            dim.interp,'linear',missing_value);
            end
         end 
         field_out(i,j,:,:)=local_field_out;
      end

   %Data is on pressure levels
   else
      fieldsize=length(size(field_out));
      switch fieldsize
         case 2
            for j=1:size(field,2)
               field_out(i,j)=interp1(dim.native.value,...
                                      squeeze(field(i,j,:)),dim.interp,...
                                      'linear',missing_value);
            end
         case 3
            if interp_dim==2
               for j=1:size(field,3)
                  field_out(i,:,j)=interp1(dim.native.value,...
                                        squeeze(field(i,:,j)),dim.interp,...
                                        'linear',missing_value);
               end
            else 
               for j=1:size(field,2)
                  field_out(i,j,:)=interp1(dim.native.value,...
                                        squeeze(field(i,j,:)),dim.interp,...
                                        'linear',missing_value);
               end
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

%Future work: this can be cleaned up with a better delineation among cases.
