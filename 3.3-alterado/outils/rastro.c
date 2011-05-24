/* Do not edit. File generated by rastro_generate. */

#include <rastro.h>

void rst_event_iiss_ptr(rst_buffer_t *ptr, u_int16_t type, u_int32_t i0, u_int32_t i1, char *s0, char *s1)
{
	rst_startevent(ptr, type<<18|0x27711);
	RST_PUT(ptr, u_int32_t, i0);
	RST_PUT(ptr, u_int32_t, i1);
	RST_PUT_STR(ptr, s0);
	RST_PUT_STR(ptr, s1);
	rst_endevent(ptr);
}

void rst_event_ii_ptr(rst_buffer_t *ptr, u_int16_t type, u_int32_t i0, u_int32_t i1)
{
	rst_startevent(ptr, type<<18|0x27700);
	RST_PUT(ptr, u_int32_t, i0);
	RST_PUT(ptr, u_int32_t, i1);
	rst_endevent(ptr);
}

void rst_event_iiiss_ptr(rst_buffer_t *ptr, u_int16_t type, u_int32_t i0, u_int32_t i1, u_int32_t i2, char *s0, char *s1)
{
	rst_startevent(ptr, type<<18|0x7771);
	RST_PUT(ptr, u_int32_t, 0x21000000);
	RST_PUT(ptr, u_int32_t, i0);
	RST_PUT(ptr, u_int32_t, i1);
	RST_PUT(ptr, u_int32_t, i2);
	RST_PUT_STR(ptr, s0);
	RST_PUT_STR(ptr, s1);
	rst_endevent(ptr);
}

