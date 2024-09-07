use libc::{c_int, c_uint};

#[repr(C)]
pub struct bsequence_s {
    // Opaque struct, implementation details hidden
}

pub type bsequence = *mut bsequence_s;

#[no_mangle]
pub extern "C" fn bsequence_create(n: c_uint) -> bsequence {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_destroy(q: bsequence) -> c_int {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_reset(q: bsequence) -> c_int {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_print(q: bsequence) -> c_int {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_push(q: bsequence, bit: c_int) -> c_int {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_circshift(q: bsequence) -> c_int {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_correlate(q: bsequence, p: bsequence) -> c_int {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_add(q: bsequence, p: bsequence) -> c_int {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_mul(q: bsequence, p: bsequence) -> c_int {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_accumulate(q: bsequence) -> c_uint {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_get_length(q: bsequence) -> c_uint {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_index(q: bsequence, i: c_uint) -> c_uint {
    unimplemented!()
}

#[no_mangle]
pub extern "C" fn bsequence_set(q: bsequence, i: c_uint, bit: c_int) -> c_int {
    unimplemented!()
}