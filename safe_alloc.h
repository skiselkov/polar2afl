/*
 * CDDL HEADER START
 *
 * This file and its contents are supplied under the terms of the
 * Common Development and Distribution License ("CDDL"), version 1.0.
 * You may only use this file in accordance with the terms of version
 * 1.0 of the CDDL.
 *
 * A full copy of the text of the CDDL should have accompanied this
 * source.  A copy of the CDDL is also available via the Internet at
 * http://www.illumos.org/license/CDDL.
 *
 * CDDL HEADER END
 */
/*
 * Copyright 2019 Saso Kiselkov. All rights reserved.
 */

#ifndef	_SAFE_ALLOC_H_
#define	_SAFE_ALLOC_H_

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifdef	__cplusplus
extern "C" {
#endif

static inline void *
safe_malloc(size_t size)
{
	void *p = malloc(size);
	if (size > 0)
		assert(p != NULL);
	return (p);
}

static inline void *
safe_calloc(size_t nmemb, size_t size)
{
	void *p = calloc(nmemb, size);
	if (nmemb > 0 && size > 0)
		assert(p != NULL);
	return (p);
}

static inline void *
safe_realloc(void *oldptr, size_t size)
{
	void *p = realloc(oldptr, size);
	if (size > 0)
		assert(p != NULL);
	return (p);
}

static inline char *
safe_strdup(const char *str2)
{
	char *str = strdup(str2);
	if (str2 != NULL)
		assert(str != NULL);
	return (str);
}

#ifdef	__cplusplus
}
#endif

#endif	/* _SAFE_ALLOC_H_ */
