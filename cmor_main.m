clear all
close all
clc

ncar_mps = parallel.importProfile('/glade/u/apps/opt/matlab/parallel/ncar_mps.mlsettings');
version=['v20230207'];

%For testing, sets to either monthly/daily file spec
file_spec_number=2;

%Script to CMORize CESM output
%Currently only functioning for the atmosphere component
%Assumes lat/lon grid, single variable timeseries (convert SE/MPAS to gridded output first)
%nadavis@ucar.edu

%Load user-supplied specification details 
cmor_specification_file = 'cmor_specification.json'; 
cmor_specification=json_load(cmor_specification_file);
 
%Load mip-specific details
[specification,specification_files]=load_specification_files([cmor_specification.cmor_specification_dir,cmor_specification.cmor_specification_format]);

%Load output file structures
output_specification_files=dir([cmor_specification.cmor_specification_dir,cmor_specification.cmor_specification_format,'*.json']);
output_specification_files=remove_entries(output_specification_files,specification_files);

%Load CESM->MIP variable dictionary
%Includes information for translating MIP requests for frequency, realm, etc.
cesm_dictionary=json_load(cmor_specification.cmor_variable_dictionary);
cesm_globals=json_load(['cesm_global_attributes.json']);
cesm_globals_names=fieldnames(cesm_globals);

%Output file info
output_specification_file=[output_specification_files(file_spec_number).folder,'/',output_specification_files(file_spec_number).name];
output_file=json_load(output_specification_file);
output=strrep(output_specification_files(file_spec_number).name,cmor_specification.cmor_specification_format,'');
output=strrep(output,'.json','');
vars=fieldnames(output_file.variable_entry);
vars_info=struct2cell(output_file.variable_entry); 
output_var_details=struct2cell(output_file.variable_entry); 
dir_input_main=[cmor_specification.case_output_dir,cmor_specification.case_name,'/'];
if strcmp(output,'Z')
   outdir='gnz';
else
   outdir='gn';
end

%Combine into single structure to pass to worker
cmor_structure=struct(...
'cmor_specification_file',cmor_specification_file,...
'cmor_specification',cmor_specification,...
'specification',specification,...
'output_specification_files',output_specification_files,...
'cesm_dictionary',cesm_dictionary,...
'cesm_globals',cesm_globals,...
'cesm_globals_names',{cesm_globals_names},...
'file_spec_number',file_spec_number,...
'dir_input_main',dir_input_main,...
'vars',{vars},...
'vars_info',{vars_info},...
'output',output,...
'output_var_details',{output_var_details},...
'output_file',output_file,...
'outdir',outdir,...
'version',version);

%Variables for job management
count=0;                         %Currently running jobs
outer_lim=-1;                    %While loop constraint; bootstrap to allow loop to begin
varsleft=length(vars);           %How many jobs left to do
varstotal=varsleft;              %How many jobs in total
varnum_current=[];               %Current job list
count_limit=2;                  %Maximum number of jobs

%Main loop
while count>outer_lim

	%Monitor each job, resubmit failures/submit new job
	while length(varnum_current) == count_limit && count>0
           vars_to_remove=[];

	   for c=1:length(varnum_current)
	      varnum=varnum_current(c);
              state=jobs{varnum}.State;
	      if strcmp(state,'finished')
	      
	         %Move on to new variable (in main loop) if successfully processed
	         count=count-1;
                 vars_to_remove=cat(1,vars_to_remove,varnum);

	      elseif strcmp(state,'failed')
	      
		      %Resubmit the same variable if the job failed
		      [jobs{varnum},vars{varnum}]=submit_batch_job(varnum,cmor_structure,cmor_specification.case_output_dir,cmor_specification.cmor_output_dir);
 	      end
	   end

           %Remove variables if successfully processed
	   if ~isempty(vars_to_remove)
	      for i=1:length(vars_to_remove)
	         varnum_current=remove_var_from_jobs(varnum_current,vars_to_remove(i));
	      end
	   end

	   %Adjust maximum number of jobs once there are no new vars to process
	   if varsleft==0
	      count_limit=length(varnum_current);
	   end
	end

	%Only submit a new variable if theres a var left to process 
	if varsleft>0

	   %Unique identifier for each job
	   varnum=get_var_num(varsleft,varstotal);
		
           [jobs{varnum},vars{varnum}]=submit_batch_job(varnum,cmor_structure,cmor_specification.case_output_dir,cmor_specification.cmor_output_dir);
	   varnum_current=add_var_to_jobs(varnum_current,varnum);
	
	   count=count+1;
	   varsleft=varsleft-1;
 
	end

        outer_lim=0;
	count=length(varnum_current);
end

%Clear out the temporary directory
cd('output')
delete('*.mat') 
delete('*.txt')
rmdir('Job*','s')

exit;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to generalize submitting a job to process a given variable
function [j,varnum]=submit_batch_job(varnum,file_structure,input_dir,output_dir)

%Set up worker information
worker_input=struct;
worker_input.varnum=varnum;
file_structure.local_var_spec=file_structure.vars_info{varnum};

%Realm, frequency information
realm=translate_cesm(file_structure.cesm_dictionary,file_structure.local_var_spec.modeling_realm,'Realms',0);
frequency=translate_cesm(file_structure.cesm_dictionary,file_structure.local_var_spec.frequency,'Frequency',0);
variables_list=translate_cesm(file_structure.cesm_dictionary,file_structure.vars{varnum},'Variables',1);
worker_input.realm=realm;
worker_input.frequency=frequency;
worker_input.variables_list=variables_list;

%Model output location
dir_input=[file_structure.dir_input_main,realm,'/proc/tseries/',frequency,'_1/'];
worker_input.dir_input=[file_structure.dir_input_main,realm,'/proc/tseries/',frequency,'_1/',file_structure.cmor_specification.case_name];

% Start PBS cluster
% Add cluster profile if not already present
if ~any(strcmp(parallel.clusterProfiles, 'ncar_mps'))
    ncar_mps = parallel.importProfile('/glade/u/apps/opt/matlab/parallel/ncar_mps.mlsettings');
end

% Start PBS cluster and submit job with custom number of workers
c = parcluster('ncar_mps');

% Setup cluster attributes - 1 node (35 workers) per cluster
jNodes = 1;
jTasks = getenv('MPSTASKS');
jWorkers = jNodes * str2num(jTasks)-1;
jAccount = getenv('MPSACCOUNT');
jQueue = getenv('MPSQUEUE');
jWalltime = getenv('MPSWALLTIME');

% Initiate job
c.ClusterMatlabRoot = getenv('NCAR_ROOT_MATLAB');
c.ResourceTemplate = append('-l select=',num2str(jNodes),':ncpus=',num2str(jTasks),':mpiprocs=', num2str(jTasks),':mem=109GB');
c.SubmitArguments = append('-A ', jAccount, ' -q ', jQueue, ' -l walltime=', jWalltime);
c.JobStorageLocation = append(getenv('PWD'), '/job_output');

% Output cluster settings
c

% Submit batch job, ensure worker has access to file so it doesn't have to create a local copy
current_dir=pwd;
function_dir=[current_dir,'/cmor_functions/'];
j = batch(c, @cmor_worker, 0, {file_structure,worker_input}, 'pool', jWorkers, 'AdditionalPaths', {dir_input,input_dir,output_dir,current_dir,function_dir});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function to get the next variable number 
function varnum=get_var_num(varsleft,varstotal)
	varnum=varstotal-varsleft+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add varnum to currently running jobs
function varnum_current=add_var_to_jobs(varnum_current,varnum)
	varnum_current=cat(1,varnum_current,varnum);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove varnum from currently running jobs
function varnum_current=remove_var_from_jobs(varnum_current,varnum)
	ind=find(varnum_current==varnum,1,'first');
	varnum_current(ind)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Remove entries from a list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function list=remove_entries(list,entries)

delete_list=[];
for i=1:length(entries)
   for j=1:length(list)
      if strcmp([list(j).folder,'/',list(j).name],entries{i})
         delete_list=cat(1,delete_list,j);
      end
   end
end
list(delete_list)=[];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load file and decode json
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decoded=json_load(file)
fid = fopen(file);
raw = fread(fid,inf);
str = char(raw');
fclose(fid);
decoded=jsondecode(str);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Load specification files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [specification,specification_files]=load_specification_files(file_preamble)

specification=struct;
specs={'CV';'fx';'grids';'coordinate';'formula_terms'};
for i=1:length(specs)
   specification_files{i}=[file_preamble,specs{i},'.json'];
   eval(['specification.',specs{i},'=json_load(''',specification_files{i},''');']);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scan dictionary and return field
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=translate_cesm(dictionary,in,field,closest_match)

out=[];

dictionary_entries=fieldnames(eval(['dictionary.',field]));
for i=1:length(dictionary_entries)
   if strcmp(dictionary_entries{i},in)
      out=eval(['dictionary.',field,'.',dictionary_entries{i}]);
   end
end

if closest_match==1
   if isempty(out)
      for i=1:length(dictionary_entries)
         if contains(in,dictionary_entries{i})
            out=eval(['dictionary.',field,'.',dictionary_entries{i}]);
         end
      end
   end
end

if isempty(out)
   out='no_match';
end

end
