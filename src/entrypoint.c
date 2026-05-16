// We need to forward routine registration from C to Rust
// to avoid the linker removing the static library.

void R_init_flod_extendr(void *dll);
void register_extendr_panic_hook(void);

void R_init_flod(void *dll) {
    register_extendr_panic_hook();
    R_init_flod_extendr(dll);
}
