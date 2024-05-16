function [ amrID ] = amr_load( filename )
% load the Chombo AMR file filename into memory
load_libamrfile();
amrID = libpointer('int32Ptr',-1);
status = libpointer('int32Ptr',-1);
calllib('libamrfile','amr_read_file',status,amrID,filename);
amr_error(status);
end

