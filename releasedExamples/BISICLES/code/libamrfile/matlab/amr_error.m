function amr_error( status )
%print an error message depending on status
    if (status.value ~= 0)
        error('amrfile:error','status = %i ',status.value);
    end
end

