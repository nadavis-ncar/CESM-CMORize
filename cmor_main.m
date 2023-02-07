clear all
close all
clc

%For testing, sets to either monthly/daily file spec
file_spec_number=4;

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
output_specification_file=[output_specification_files(i).folder,'/',output_specification_files(i).name];
output_file=json_load(output_specification_file);
vars=fieldnames(output_file.variable_entry);

%Combine into single structure to pass to worker
cmor_structure=struct(...
'cmor_specification_file',cmor_specification_file,...
'cmor_specification',cmor_specification,...
'specification',specification,...
'output_specification_files',output_specification_files,...
'cesm_dictionary',cesm_dictionary,...
'cesm_globals',cesm_globals,...
'cesm_globals_names',cesm_globals_names,...
'file_spec_number',file_spec_number);

%Variables for job management
count=0;                         %Currently running jobs
outer_lim=-1;                    %While loop constraint; bootstrap to allow loop to begin
varsleft=length(vars);           %How many jobs left to do
varstotal=varsleft;              %How many jobs in total
varnum_current=[];               %Current job list
count_limit=4;                  %Maximum number of jobs

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
		      [jobs{varnum},files{varnum}]=submit_batch_job(varnum,cmor_structure,cmor_specification.case_input_dir,cmor_specification.cmor_output_dir);
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
		
      [jobs{varnum},files{varnum}]=submit_batch_job(varnum,cmor_structure,cmor_specification.case_input_dir,cmor_specification.cmor_output_dir);
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
function [j,input_file]=submit_batch_job(varnum,file_structure,input_dir,output_dir)

%Set up worker information
worker_input=struct;
worker_input.varnum=varnum;

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
%c.NumWorkers=str2num(jTasks);

% Initiate job
c.ClusterMatlabRoot = getenv('NCAR_ROOT_MATLAB');
c.ResourceTemplate = append('-l select=',num2str(jNodes),':ncpus=',num2str(jWorkers),':mpiprocs=', num2str(jWorkers),':mem=109GB');
c.SubmitArguments = append('-A ', jAccount, ' -q ', jQueue, ' -l walltime=', jWalltime);
c.JobStorageLocation = append(getenv('PWD'), '/output');

% Output cluster settings
c

% Submit batch job, ensure worker has access to file so it doesn't have to create a local copy
current_dir=pwd;
j = batch(c, @cmor_worker, 0, {file_structure,worker_input}, 'pool', jWorkers, 'AdditionalPaths', {input_dir;output_dir;current_dir});

input_file=worker_input.file

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
