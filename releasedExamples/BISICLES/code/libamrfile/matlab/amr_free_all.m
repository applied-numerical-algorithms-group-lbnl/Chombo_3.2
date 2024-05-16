function amr_free_all( )
%free all memory holding Chombo AMR data
load_libamrfile();
calllib('libamrfile','amr_free_all');
end
