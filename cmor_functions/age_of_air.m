%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate age of air
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tracer=age_of_air(tracer,p,y,p_ref,y_ref)

tracer=movmean(tracer,12,4);

%Reference time series
[val,y_ref]=min(abs(y-y_ref));
[val,p_ref]=min(abs(p-p_ref));
reference_timeseries=squeeze(mean(tracer(:,y_ref,p_ref,:),1));

%Calculation
parfor i=1:size(tracer,1)
   tracer_local=squeeze(tracer(i,:,:,:));
   aoa_out_local=zeros(size(tracer_local));
   reference_timeseries_local=repmat(reference_timeseries(:),[1 length(reference_timeseries)]);
   for j=1:size(tracer_local,1)
      for k=1:size(tracer_local,2)
         timeseries=squeeze(tracer_local(j,k,:));
         timeseries=repmat(timeseries(:),[1 length(timeseries)]);
         [val,time_index]=min(abs(timeseries-reference_timeseries_local),[],2);
         aoa_out_local(j,k,:)=time_index'-([1:size(timeseries,1)]/12);
      end
   end
   tracer(i,:,:,:)=aoa_out_local;
end
