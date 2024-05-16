function [ time ] = amr_query_time( amrID )
%lookup the time 
load_libamrfile();
status = libpointer('int32Ptr',-1);
timep=libpointer('doublePtr',0.0);
calllib('libamrfile','amr_query_time',status,timep,amrID);
amr_error(status);
time=timep.value;
end
