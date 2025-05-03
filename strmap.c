#include "msppriv.h"
#include "khashl.h"
KHASHL_MAP_INIT(KH_LOCAL, strmap_t, strmap, const char*, int32_t, kh_hash_str, kh_eq_str)

int msp_verbose = 3;

/**************
 * String map *
 **************/

msp_strmap_t *msp_strmap_init(void)
{
	msp_strmap_t *m;
	m = MSP_CALLOC(msp_strmap_t, 1);
	m->h = strmap_init();
	return m;
}

void msp_strmap_destroy(msp_strmap_t *m)
{
	int32_t i;
	strmap_destroy((strmap_t*)m->h);
	for (i = 0; i < m->n; ++i)
		free(m->a[i]);
	free(m->a);
	free(m);
}

int32_t msp_strmap_add(msp_strmap_t *m, const char *s)
{
	strmap_t *h = (strmap_t*)m->h;
	int absent;
	khint_t k;
	k = strmap_put(h, s, &absent);
	if (absent) {
		MSP_GROW(char*, m->a, m->n, m->m);
		kh_key(h, k) = m->a[m->n] = msp_strdup(s);
		kh_val(h, k) = m->n++;
	}
	return kh_val(h, k);
}

int32_t msp_strmap_get(const msp_strmap_t *m, const char *s)
{
	const strmap_t *h = (const strmap_t*)m->h;
	khint_t k;
	k = strmap_get(h, s);
	return k != kh_end(h)? kh_val(h, k) : -1;
}
