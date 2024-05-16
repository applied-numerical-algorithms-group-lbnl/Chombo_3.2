function  load_libamrfile()
%load libamrfile shared library if required
    if (~libisloaded('libamrfile'))
        loadlibrary('libamrfile.so','src/libamrfile.H');
    end
end

