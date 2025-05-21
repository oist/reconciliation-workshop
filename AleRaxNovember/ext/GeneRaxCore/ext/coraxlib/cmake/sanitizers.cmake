# Adds the compiler and linker flags necessary to enable ASan to the given target.
function (target_add_address_sanitizer TARGET)
    target_compile_options(${TARGET} PRIVATE -fsanitize=address)
    target_link_options(${TARGET} PRIVATE -fsanitize=address)
endfunction ()
