function cmor_worker(cmor_structure,worker_input)

disp('Starting')
cmor_structure
worker_input

%Loading necessary worker, cmor information from input structures
cmor_specification_file=cmor_structure.cmor_specification_file;
cmor_specification=cmor_structure.cmor_specification;
specification=cmor_structure.specification;
output_specification_files=cmor_structure.output_specification_files;
cesm_dictionary=cmor_structure.cesm_dictionary;
cesm_globals=cmor_structure.cesm_globals;
cesm_globals_names=cmor_structure.cesm_globals_names;
file_spec_number=cmor_structure.file_spec_number;
v=worker_input.varnum;
output_file=cmor_structure.output_file;
output=cmor_structure.output;
vars=cmor_structure.vars;
vars_info=cmor_structure.vars_info;
output_var_details=cmor_structure.output_var_details;
local_var_spec=cmor_structure.local_var_spec;
dir_input_main=cmor_structure.dir_input_main;
realm=worker_input.realm;
frequency=worker_input.frequency;
variables_list=worker_input.variables_list;
dir_input=worker_input.dir_input;
version=cmor_structure.version;

tem=0;

tic
disp(vars{v})

local_var_spec.averaging=averaging_flags(local_var_spec.cell_methods);

%Gather variable info, special operations, scale_factors; or just load the variable 
if ~isstruct(variables_list) 
   variables.var{1}=variables_list;
   variables.operation=[];
   variables.eval=[];
   variables.option=[];
else
   [variables]=translate_cesm_variable(variables_list);
end
load_ps=0;
save_ps=0;
load_lev=0;

%Grab global variables for the file
globals=global_attributes(specification.CV.CV.required_global_attributes,specification.CV.CV,output,cmor_specification,specification,frequency,vars{v},cesm_globals);
globals=set_index_variables(globals,cmor_specification.cmor_case_name);

%For ambiguous variables where there are two or more supplied output options
if strcmp(fieldnames(variables),'option')
   for j=1:length(variables.option)
      if contains(local_var_spec.comment,variables.option(j).string)
         variables.var{1}=variables.option(j).var;
      end
   end
end

%Set TEM variable(s)
if isempty(variables.var) & strcmp(variables.operation,'TEM') & tem==0
   disp('TEM variable requested, now calculating all TEM output variables')
   variables.var{1}='Uzm';
   variables.var{2}='Vzm';
   variables.var{3}='Wzm';
   variables.var{4}='VTHzm';
   variables.var{5}='UVzm';
   variables.var{6}='UWzm';
   variables.var{7}='THzm'; 
   tem=1;
end

%Check if the variable will be output or not
if ~isempty(variables.var) & ~strcmp(variables.var{1},'no_match')

   %Assemble input file list(s)
   for vin=1:length(variables.var)
      var_files=dir([dir_input,'*.',variables.var{vin},'.*']);
      for i=1:length(var_files)
         file_list{vin}{i}=[var_files(i).folder,'/',var_files(i).name];
      end
   end

   %Output grid label
   grid_label='gn';
   if strcmp(output(end),'Z')
      grid_label=[grid_label,'z'];
   end

   %MIP-compatible model name
   for j=1:length(globals)
      if strcmp(globals(j).name,'source_id')
         id_index=j;
      end
   end

   if strcmp(frequency,'month')
      date_ind=15;
   else
      date_ind=19;
   end

   %If there are some files done, then trim the file list
   if ~strcmp(variables.operation,'TEM')   
      dir_output=[cmor_specification.cmor_output_dir,cmor_specification.case_name,'/postprocess/output/',output,'/',vars{v},'/',grid_label,'/',version,'/'];
      if ~exist(dir_output)
         mkdir(dir_output)
      end
      files_done=[];
      for i=1:length(file_list{1})
         file_name=file_list{1}{i};
         outfile=[dir_output,vars{v},'_',output,'_',globals(id_index).value,'_',cmor_specification.cmor_experiment,...
                  '_',cmor_specification.cmor_case_name,'_',grid_label,'_',file_name(end-date_ind:end-3),'.nc'];
         if isfile(outfile)
            for vin=1:length(variables.var)
               file_list{vin}(i)=[];
            end
         end
      end
   end

   %The odd case where everything is already done
   if isempty(file_list{1})
      disp('No data left to process for this variable')
   end

   %Initiate postprocessing
   while ~isempty(file_list{1})
      %Load variable(s) 
      for vin=1:length(variables.var)
         if vin>1
           %Ensure file lists are consistent dates
           file_name=file_list{1}{1};
           date_str=file_name(end-date_ind:end-3);
           file_name=file_list{vin}{1};
           date_str_2=file_name(end-date_ind:end-3);
           f_2=1;
           while ~strcmp(date_str,date_str_2)
              f_2=f_2+1;
              if f_2>length(file_list{vin})
                 error('No files with dates matching main variable file')
              else
                 file_name=file_list{vin}{f_2};
                 date_str_2=file_name(end-date_ind:end-3);
              end
           end
           file_name=file_list{vin}{f_2};
         else
           file_name=file_list{vin}{1};
         end            

         date_string=file_name(end-date_ind:end-3);
         date=ncread(file_name,'date');
         time=ncread(file_name,'time');
         dims=parse_string(local_var_spec.dimensions);

         %Get dimensions of input variable
         %Deal with useless singular dimension "zlon"
         remove_list=[];
         file_structure=ncinfo(file_name);
         for j=1:length(file_structure.Variables)
            if strcmp(file_structure.Variables(j).Name,variables.var{vin})
               for k=1:length(ncinfo(file_name).Variables(j).Dimensions)
                  if strcmp(ncinfo(file_name).Variables(j).Dimensions(k).Name,'zlon') 
                     remove_list=k;
                  else
                     dim_native_values{k}=ncread(file_name,ncinfo(file_name).Variables(j).Dimensions(k).Name);
                     dim_native_names{k}=ncinfo(file_name).Variables(j).Dimensions(k).Name;
                     dim_native(k)=k;
                  end
               end       
            end
         end
         dim_native_names(remove_list)=[];
         dim_native(remove_list)=[];
         dim_native=dim_native-length(remove_list);
         dim_native_names_list=dim_native_names;

         %Load dimension info, connect output dims to input dims
         for j=1:length(dims)

           dim{vin}{j}.native.name=translate_cesm(cesm_dictionary,dims{j},'Dimension',1);   

            %Single-level/subdomain 
            if strcmp(dim{vin}{j}.native.name,'no_match')
               dim{vin}{j}.native.name=translate_cesm(cesm_dictionary,eval(['specification.coordinate.axis_entry.',dims{j},'.out_name']),'Dimension',1);
            end

            %Already interpolated?
            if strcmp(dim{vin}{j}.native.name,'lev')
               for k=1:length(dim_native_names)
                  if strcmp(dim_native_names{k},'ilev')
                     dim{vin}{j}.native.name='ilev';
                     load_lev=1;
                  end
               end
            end

            if ~strcmp(dim{vin}{j}.native.name,'no_match')
               dim_native(strcmp(dim{vin}{j}.native.name,dim_native_names))=[];
               dim_native_names(strcmp(dim{vin}{j}.native.name,dim_native_names))=[]; 
               dim{vin}{j}.native.value=ncread(file_name,dim{vin}{j}.native.name);
               dim{vin}{j}.out.name=dims{j};
               dim{vin}{j}.out.info=dimension_info(specification,dims{j},cesm_dictionary);
               dim{vin}{j}.interp_direction='';
               dim{vin}{j}.interp=[];

               %Assimilate MIP lev info for native output
               if strcmp(dim{vin}{j}.out.info,'lev')
                  dim{vin}{j}.out.info=dimension_info(specification,'standard_hybrid_sigma',cesm_dictionary);
               end

               %If we need to interpolate to a new vertical grid
               if (strcmp('lev',dim{vin}{j}.native.name) | strcmp('ilev',dim{vin}{j}.native.name)) ...
                  & sum(strcmp(fieldnames(specification.coordinate.axis_entry),dim{vin}{j}.out.name))>0
                  dim{vin}{j}=get_requested_axis_value(dim{vin}{j},specification.coordinate.axis_entry);
                  dim{vin}{j}.interp_direction='vertical';
               end     

               %Convert from hPa to Pa
               if strcmp(dim{vin}{j}.native.name,'ilev') | strcmp(dim{vin}{j}.native.name,'lev')
                  dim{vin}{j}.native.value=dim{vin}{j}.native.value*100;
               end

               %Set load settings if we need to interpolate the vertical grid
               if strcmp(dim{vin}{j}.native.name,'lev')                
                  lev_dim=j;
                  load_lev=1;
                  load_ps=1;
                  if isempty(dim{vin}{j}.interp)
                     save_ps=1;
                  end
                  if contains(dim{vin}{j}.out.name,'lev')              
                     dim{vin}{j}.out.name='lev';
                     dim{vin}{j}.interp='';
                  end
               end

               %Add consistent time units
               if strcmp(dim{vin}{j}.native.name,'time')
                  dim{vin}{j}.out.info.units=ncreadatt(file_name,'time','units');
               end

               %Case where the output dimension is really a single-level field
               %Could be a creative way to handle this but it seems to be a single case
            elseif strcmp(dims{j},'height2m')
               dim{vin}{j}.out.name=dims{j};
               dim{vin}{j}.out.info=dimension_info(specification,dims{j},cesm_dictionary);
               dim{vin}{j}.native.value=2.0;
               dim{vin}{j}.native.name='height2m';
               dim{vin}{j}.interp_direction='';
               dim{vin}{j}.interp=[];
               dim_native_names_list{j}='height2m';
               dim_native_values{j}=2.0;
            else
               error('unrecognized dimension name requested by MIP')
            end
         end

        %Is lev dimension being removed (integration, single-level interp)? Need to load native grid info
         if ~isempty(dim_native)
            if strcmp(dim_native_names,'lev')
               load_lev=1;
               load_ps=1;
            end 
         end

         %Load lev formula fields
         if load_lev==1  
            a=ncread(file_name,'hyam')*1e5;
            b=ncread(file_name,'hybm');
            ai=ncread(file_name,'hyai')*1e5;
            bi=ncread(file_name,'hybi');
            ilev=ncread(file_name,'ilev');
            pure_pressure=find(b==0,1,'last');
         end

         %Load variable
         var{vin}.native=struct;
         var{vin}.native.value=ncread(file_name,variables.var{vin});
         if strcmp(variables.operation,'TEM') | ~isempty(remove_list) 
            var{vin}.native.value=squeeze(var{vin}.native.value);
         end
         var{vin}.native.name=variables.var{vin};

         %Load ps from equivalent-date file
         if load_ps==1
            var_files=dir([dir_input,'*.PS.*',date_string,'*']);
            file_name_ps=[var_files(1).folder,'/',var_files(1).name];
            ps_out.value=ncread(file_name_ps,'PS');
            for j=1:length(vars)
               if strcmp(vars{j},'ps')
                  ps_out.info=vars_info{j};
               end
            end
         end 
      end

      %Inherit output dimensions of input var
      var_out.dim=dim{1};

      %Note dimensions to be removed 
      var_out.dim_names_to_remove=dim_native_names;
      var_out.dim_to_remove=dim_native;

      %Perform arithmetic operations, or load secondary variable for a calculation
      %Eval converts a string from JSON to code and executes
      if isstruct(variables)
         if ~isempty(variables.eval) 
            eval(['var_out.native.value=',variables.eval{1},';'])
         else
            var_out.native=var{1}.native;
            if length(var)>1
               var_out.native_secondary=var{2}.native;
            end
         end 
      else
         var_out.native=var{1}.native;
      end

      if strcmp(variables.operation,'TEM')
         for j=1:length(variables.var)
            var_out.var{j}=var{j}.native.value;
         end
      end

      clear var

      for j=1:length(dim_native_names_list)
         var_out.native_dim(j).name=dim_native_names_list{j};
         var_out.native_dim(j).value=dim_native_values{j};
      end

      %Special operations
      if ~isempty(variables.operation)
         operation=variables.operation;
         switch operation
            case 'age_of_air'
               var_out.native.value=age_of_air(var_out.native.value,var_out.dim{3}.native.value,var_out.dim{2}.native.value,variables.p0,variables.lat0);
            case 'sum'
               var_out=calculate_sum(var_out,variables.axis);
            case 'ozone_integral'
               var_out.native.value=convert_units(var_out.native.value,ps_out.value,a,b,'molar_mixing_ratio_to_DU');
               if strcmp(fieldnames(var_out),'native_secondary')
                  var_out.native.value=calculate_integral(var_out.native.value,ps_out.value,a,b,var_out.native_secondary);
               else
                  var_out.native.value=calculate_integral(var_out.native.value,ps_out.value,a,b);
               end
            case 'burden'
               var_out.native.value=var_out.native.value.*var_out.native_secondary.value;
disp('burden')
               var_out.native.value=convert_units(var_out.native.value,ps_out.value,a,b,'number_density_to_molar_mixing_ratio');
disp('integral')
               var_out.native.value=calculate_integral(var_out.native.value,ps_out.value,a,b);
disp('done')
            case 'max_value'
               var_out=calculate_maximum(var_out,variables.axis);
            case 'omega_to_w'
               var_out.native.value=var_out.native.value./var_out.native_secondary.value;
               var_out.native.value=omega_to_w(var_out.native.value,ps_out.value,a,b);
            case 'TEM'
               var_out=calculate_tem(var_out,str2num(output_file.Header.missing_value));
         end
      end

      %Ensure dimension ordering is correct, if not, map out the reorder
      permute_order=[];
      for j=1:length(var_out.dim)
         if strcmp(var_out.dim{j}.native.name,dim_native_names_list{j})
            permute_order=[permute_order j];
         else
            for k=1:length(dim_native_names_list)
               if strcmp(var_out.dim{j}.native.name,dim_native_names_list{k})
                  permute_order=[permute_order k];
               end
            end
         end
      end
disp('reordering')
      %Only attempt reordering if necessary - typically for single pressure-level output
      if length(permute_order) == length(size(var_out.native.value))    
         if sum(diff(permute_order)<0)>0
            if max(permute_order) > length(size(var_out.native.value))
               [val,loc]=max(permute_order);
               permute_order(loc)=permute_order(loc)-1;
            end
            var_out.native.value=permute(var_out.native.value,permute_order);
            for k=1:length(var_out.native_dim)
               native_names_reorder{k}=var_out.native_dim(permute_order(k)).name;
            end
            for k=1:length(var_out.native_dim)
               var_out.native_dim(k).name=native_names_reorder{k};
            end
         end
      end
disp('interp')
      %Do any interpolation
      for j=1:length(var_out.dim) 
         if ~isempty(var_out.dim{j}.interp)
            %Assign dim index to be interpolated - mapping between native and output space
            if strcmp(var_out.dim{j}.interp_direction,'vertical')
               dimsearch='lev';
            else
               dimsearch=var_out.dim{j}.name;
            end
            interp_dim=get_dim_index(var_out.native_dim,dimsearch,'contains');
            %Carry out interpolation
            if strcmp(var_out.dim{j}.interp_direction,'vertical') 
               if strcmp(variables.operation,'TEM')
                  tem_output_vars=fieldnames(var_out.tem);
                  for k=1:length(tem_output_vars)
                     eval(['var_out.tem.',tem_output_vars{k},'=interpolate_field(var_out.tem.',tem_output_vars{k},',interp_dim,var_out.dim{j},NaN);']);
                  end
               elseif load_ps==1
                  var_out.native.value=interpolate_field(var_out.native.value,interp_dim,var_out.dim{j},NaN,a,b,ps_out.value);
               else
                  var_out.native.value=interpolate_field(var_out.native.value,interp_dim,var_out.dim{j},NaN);
               end
            else
               var_out.native.value=interpolate_field(var_out.native.value,interp_dim,var_out.dim{j},NaN);
            end
            var_out.dim{j}.out.value=var_out.dim{j}.interp;
         else
            var_out.dim{j}.out.value=var_out.dim{j}.native.value;
         end
      end

      %Averaging - could generalize to any dim (will this ever need to change?)
      for j=1:length(local_var_spec.averaging)
         if ~strcmp(local_var_spec.averaging{j},'time') & strcmp(local_var_spec.averaging{j},'longitude') & length(remove_list)==0
            var_out.native.value=squeeze(mean(var_out.native.value,1,'omitnan'));
            var_out.native.value(isnan(var_out.native.value))=str2num(output_file.Header.missing_value);
            ps_out.value=mean(ps_out.value,1);
         end
      end

      %Ensure dimension directions are consistent with file specifications
      var_out=stored_direction(var_out);   

      do_ab=0;
      %Bnd variables
      for j=1:length(var_out.dim)
         if strcmp(var_out.dim{j}.out.name,'lev')
            var_out.dim{j}.out.bnds = create_bnds(var_out.dim{j}.out,ilev);
            do_ab=j;
         elseif strcmp(var_out.dim{j}.out.name,'time')
            var_out.dim{j}.out.bnds= ncread(file_name,'time_bnds');
         elseif length(var_out.dim{j}.out.value)>1
            var_out.dim{j}.out.bnds = create_bnds(var_out.dim{j}.out);
         end
      end

      %Bnd variables (vertical)
      if do_ab~=0
         var_out.dim{do_ab}.out.formula.a_bnds.value = create_bnds(a,ai);
         var_out.dim{do_ab}.out.formula.a_bnds.info = get_formula_variable(specification.formula_terms.formula_entry,'a_bnds');
         var_out.dim{do_ab}.out.formula.b_bnds.value = create_bnds(b,bi);
         var_out.dim{do_ab}.out.formula.b_bnds.info = get_formula_variable(specification.formula_terms.formula_entry,'b_bnds');
         var_out.dim{do_ab}.out.formula.a.value=a;
         var_out.dim{do_ab}.out.formula.a.info = get_formula_variable(specification.formula_terms.formula_entry,'a');
         var_out.dim{do_ab}.out.formula.b.value=b;
         var_out.dim{do_ab}.out.formula.b.info = get_formula_variable(specification.formula_terms.formula_entry,'b');
         if var_out.dim{do_ab}.flip==1
            var_out.dim{do_ab}.out.formula.a_bnds.value=flip(var_out.dim{do_ab}.out.formula.a_bnds.value,2);
            var_out.dim{do_ab}.out.formula.b_bnds.value=flip(var_out.dim{do_ab}.out.formula.b_bnds.value,2);
            var_out.dim{do_ab}.out.formula.a.value=flip(var_out.dim{do_ab}.out.formula.a.value,1);
            var_out.dim{do_ab}.out.formula.b.value=flip(var_out.dim{do_ab}.out.formula.b.value,1);
         end
      end

      %Set info for NetCDF
      var_out.info=output_var_details{v};
disp('save')
      %Save to file      
      if strcmp(variables.operation,'TEM')
         tem_vars=struct2cell(var_out.tem);
         for k=1:length(fieldnames(var_out.tem))
            var_out.native.value=eval(['var_out.tem.',tem_output_vars{k},';']);
            var_out.info=output_var_details{find(strcmp(vars,tem_output_vars{k}),1,'first')};
            dir_output=[cmor_specification.cmor_output_dir,cmor_specification.case_name,'/postprocess/output/',output,'/',tem_output_vars{k},'/',grid_label,'/',version,'/'];
            if ~exist(dir_output)
               mkdir(dir_output)
            end
            var_out=set_fill_values(var_out,str2num(output_file.Header.missing_value));
            outfile=[dir_output,tem_output_vars{k},'_',output,'_',globals(id_index).value,'_',cmor_specification.cmor_experiment,...
                    '_',cmor_specification.cmor_case_name,'_',grid_label,'_',file_name(end-date_ind:end-3),'.nc'];
            create_netcdf(var_out,outfile,globals);
         end
      else
         var_out=set_fill_values(var_out,str2num(output_file.Header.missing_value));
         outfile=[dir_output,vars{v},'_',output,'_',globals(id_index).value,'_',cmor_specification.cmor_experiment,...
                  '_',cmor_specification.cmor_case_name,'_',grid_label,'_',file_name(end-date_ind:end-3),'.nc'];
         if save_ps==1
            create_netcdf(var_out,outfile,globals,ps_out);
         else
            create_netcdf(var_out,outfile,globals);
         end
      end 

      %Remove files from list 
      for k=1:length(variables.var)
         file_list{k}(1)=[];
      end

      clear dim dim_native var_out ps_out dim_out
     
      toc
   end
end
