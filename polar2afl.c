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

#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "avl.h"
#include "safe_alloc.h"

#ifndef	MIN
#define	MIN(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef	MAX
#define	MAX(x, y)	((x) > (y) ? (x) : (y))
#endif

typedef struct {
	double		alpha;
	double		Cl;
	double		Cd;
	double		Cm;
	avl_node_t	node;
} polar_t;

typedef struct {
	double		Re;
	char		*afl_data;
	avl_tree_t	polars;
	avl_node_t	node;
} polar_diag_t;

typedef struct {
	char		preamble[1024];
	avl_tree_t	diags;
} afl_t;

static double	blend_range = 10;	/* degrees */

static void afl_free(afl_t *afl);

static double
wavg(double x, double y, double w)
{
	assert(w >= 0.0 && w <= 1.0);
	return (x + (y - x) * w);
}

static inline double
clamp(double x, double min_val, double max_val)
{
	if (min_val < max_val) {
		if (x < min_val)
			return (min_val);
		if (x > max_val)
			return (max_val);
	} else {
		if (x > min_val)
			return (min_val);
		if (x < max_val)
			return (max_val);
	}
	return (x);
}

static inline double
iter_fract(double x, double min_val, double max_val, bool clamp_output)
{
	x = (x - min_val) / (max_val - min_val);
	if (clamp_output)
		x = clamp(x, 0, 1);
	return (x);
}

/*
 * Removes all leading & trailing whitespace from a line.
 */
void
strip_space(char *line)
{
	char	*p;
	size_t	len = strlen(line);

	/* strip leading whitespace */
	for (p = line; *p != 0 && isspace(*p); p++)
		;
	if (p != line) {
		memmove(line, p, (len + 1) - (p - line));
		len -= (p - line);
	}

	for (p = line + len - 1; p >= line && isspace(*p); p--)
		;
	p[1] = 0;
}

static int
polar_diag_compar(const void *a, const void *b)
{
	const polar_diag_t *pa = a, *pb = b;

	if (pa->Re < pb->Re)
		return (-1);
	if (pa->Re > pb->Re)
		return (1);
	return (0);
}

static int
polar_compar(const void *a, const void *b)
{
	const polar_t *pa = a, *pb = b;

	if (pa->alpha < pb->alpha)
		return (-1);
	if (pa->alpha > pb->alpha)
		return (1);
	return (0);
}

static void
polar_diag_free(polar_diag_t *diag)
{
	void *cookie = NULL;
	polar_t *polar;

	free(diag->afl_data);
	while ((polar = avl_destroy_nodes(&diag->polars, &cookie)) != NULL)
		free(polar);
	avl_destroy(&diag->polars);
	free(diag);
}

static polar_diag_t *
polar_diag_alloc(void)
{
	polar_diag_t *diag = safe_calloc(1, sizeof (*diag));
	avl_create(&diag->polars, polar_compar, sizeof (polar_t),
	    offsetof(polar_t, node));
	return (diag);
}

static polar_diag_t *
afl_parse_polar_diag(FILE *fp, const char *filename, int *line_nr)
{
	char *line = NULL;
	size_t linecap = 0;
	int n;
	polar_diag_t *diag = polar_diag_alloc();

	n = getline(&line, &linecap, fp);
	/* Probably an EOF, just exit */
	if (n <= 0)
		goto errout;
	if (sscanf(line, "%lf", &diag->Re) != 1) {
		fprintf(stderr, "%s:%d: error parsing Reynolds number\n",
		    filename, *line_nr);
		goto errout;
	}
	diag->Re *= 1000000.0;
	(*line_nr)++;
	diag->afl_data = safe_strdup(line);
	/* skip the "alpha cl cd cm:" line */
	if (getline(&line, &linecap, fp) <= 0) {
		fprintf(stderr, "%s:%d: premature EOF\n", filename, *line_nr);
		goto errout;
	}
	(*line_nr)++;

	for (;getline(&line, &linecap, fp) > 0; (*line_nr)++) {
		double alpha, Cl, Cd, Cm;
		polar_t *polar;
		avl_index_t where;

		if (sscanf(line, "%lf %lf %lf %lf", &alpha, &Cl, &Cd, &Cm) !=
		    4) {
			fprintf(stderr, "%s:%d: error parsing polar line\n",
			    filename, *line_nr);
			goto errout;
		}
		polar = safe_calloc(1, sizeof (*polar));
		polar->alpha = alpha;
		polar->Cl = Cl;
		polar->Cd = Cd;
		polar->Cm = Cm;
		if (avl_find(&diag->polars, polar, &where) != NULL) {
			fprintf(stderr, "%s:%d: duplicate alpha detected\n",
			    filename, *line_nr);
			goto errout;
		}
		avl_insert(&diag->polars, polar, where);
		/* The 180 degree polar is the last one in a diagram */
		if (alpha >= 180.0)
			break;
	}

	free(line);
	return (diag);
errout:
	polar_diag_free(diag);
	free(line);
	return (NULL);
}

static afl_t *
afl_alloc(void)
{
	afl_t *afl = safe_calloc(1, sizeof (*afl));

	avl_create(&afl->diags, polar_diag_compar,
	    sizeof (polar_diag_t), offsetof(polar_diag_t, node));
	return (afl);
}

static afl_t *
afl_parse(const char *filename)
{
	afl_t *afl;
	char *line = NULL;
	size_t linecap = 0;
	FILE *fp;
	int line_nr;
	polar_diag_t *diag;

	fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "Can't open %s: %s\n", filename,
		    strerror(errno));
		return (NULL);
	}

	afl = afl_alloc();

	for (line_nr = 0; getline(&line, &linecap, fp) > 0; line_nr++) {
		if (line_nr < 19) {
			strncat(afl->preamble, line,
			    sizeof (afl->preamble) - 1);
			continue;
		}
		/*
		 This is just the number of diags in the file, we don't
		 * need that, we parse until we can't get any more polars
		 */
		if (line_nr == 19) {
			line_nr++;
			break;
		}
	}
	while ((diag = afl_parse_polar_diag(fp, filename, &line_nr)) != NULL) {
		avl_index_t where;

		if (avl_find(&afl->diags, diag, &where) != NULL) {
			fprintf(stderr, "%s: duplicate polar diagram for "
			    "Re = %f\n", filename, diag->Re);
			polar_diag_free(diag);
			goto errout;
		}
		avl_insert(&afl->diags, diag, where);
	}

	free(line);
	fclose(fp);
	return (afl);
errout:
	afl_free(afl);
	free(line);
	fclose(fp);
	return (NULL);
}

static void
afl_write_diag(FILE *fp, const polar_diag_t *diag)
{
	if (avl_numnodes(&diag->polars) == 0)
		return;

	if (diag->afl_data != NULL)
		fputs(diag->afl_data, fp);
	fprintf(fp, "alpha cl cd cm:\n");

	for (double alpha = -180; alpha <= 180;) {
		double step = ((alpha < -20.0 || alpha >= 20.0) ? 1.0 : 0.1);
		polar_t srch = { .alpha = alpha };
		polar_t *p1, *p2;
		avl_index_t where;
		double Cl, Cd, Cm;

		p1 = avl_find(&diag->polars, &srch, &where);
		if (p1 != NULL) {
			Cl = p1->Cl;
			Cd = p1->Cd;
			Cm = p1->Cm;
		} else {
			p1 = avl_nearest(&diag->polars, where, AVL_BEFORE);
			p2 = avl_nearest(&diag->polars, where, AVL_AFTER);

			assert(p1 != NULL || p2 != NULL);
			if (p1 != NULL && p2 != NULL) {
				double w = iter_fract(alpha, p1->alpha,
				    p2->alpha, true);

				Cl = wavg(p1->Cl, p2->Cl, w);
				Cd = wavg(p1->Cd, p2->Cd, w);
				Cm = wavg(p1->Cm, p2->Cm, w);
			} else {
				if (p1 == NULL)
					p1 = p2;
				Cl = p1->Cl;
				Cd = p1->Cd;
				Cm = p1->Cm;
			}
		}

		fprintf(fp, "%6.1f %8.5f %8.5f %8.5f\n", alpha, Cl, Cd, Cm);

		alpha = round((alpha + step) * 10.0) / 10.0;
	}


}

static bool
afl_write(const afl_t *afl, const char *filename)
{
	FILE *fp = fopen(filename, "w");

	if (fp == NULL) {
		fprintf(stderr, "Can't write airfoil %s: %s\n", filename,
		    strerror(errno));
		return (false);
	}

	if (afl->preamble != NULL)
		fprintf(fp, "%s", afl->preamble);
	fprintf(fp, "%d\n", (int)avl_numnodes(&afl->diags));

	for (const polar_diag_t *diag = avl_first(&afl->diags);
	    diag != NULL; diag = AVL_NEXT(&afl->diags, diag)) {
		afl_write_diag(fp, diag);
	}

	fclose(fp);
	return (true);
}

static void
afl_free(afl_t *afl)
{
	void *cookie = NULL;
	polar_diag_t *diag;

	while ((diag = avl_destroy_nodes(&afl->diags, &cookie)) !=
	    NULL) {
		polar_diag_free(diag);
	}
	avl_destroy(&afl->diags);
	free(afl);
}

static bool
xfoil_polar_parse(afl_t *xfoil, const char *filename)
{
	FILE *fp;
	polar_diag_t *diag;
	bool header = true;
	char *line = NULL;
	size_t linecap = 0;
	avl_index_t where;

	fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "Error opening %s: %s\n", filename,
		    strerror(errno));
		return (false);
	}
	diag = polar_diag_alloc();

	while (getline(&line, &linecap, fp) > 0) {
		strip_space(line);
		if (strlen(line) == 0)
			continue;
		if (header) {
			const char *p;
			double mantissa, exponent;

			if ((p = strstr(line, "Re =")) != NULL) {
				if (sscanf(&p[4], "%lf e %lf", &mantissa,
				    &exponent) != 2) {
					fprintf(stderr, "%s: error parsing "
					    "Re\n", filename);
					goto errout;
				}
				diag->Re = mantissa * pow(10, exponent);
			} else if ((p = strstr(line, "-------")) != NULL) {
				header = false;
			}
		} else {
			polar_t *polar = safe_calloc(1, sizeof (*polar));
			avl_index_t where;
			double CDp;

			if (sscanf(line, "%lf %lf %lf %lf %lf",
			    &polar->alpha, &polar->Cl, &polar->Cd, &CDp,
			    &polar->Cm) != 5) {
				fprintf(stderr, "%s: error parsing polar\n",
				    filename);
				goto errout;
			}
			if (avl_find(&diag->polars, polar, &where) != NULL) {
				free(polar);
				continue;
			}
			avl_insert(&diag->polars, polar, where);
		}
	}
	if (avl_find(&xfoil->diags, diag, &where) != NULL) {
		fprintf(stderr, "%s: is a duplicate of an existing polar "
		    "at Re = %fM\n", filename, diag->Re / 1000000);
		goto errout;
	}
	avl_insert(&xfoil->diags, diag, where);

	free(line);
	fclose(fp);
	return (true);
errout:
	polar_diag_free(diag);
	free(line);
	fclose(fp);
	return (false);
}

static polar_t *
find_nearest_polar(const polar_diag_t *diag, double alpha)
{
	polar_t srch = { .alpha = alpha };
	avl_index_t where;
	polar_t *polar = avl_find(&diag->polars, &srch, &where);

	if (polar == NULL)
		polar = avl_nearest(&diag->polars, where, AVL_BEFORE);
	if (polar == NULL)
		polar = avl_nearest(&diag->polars, where, AVL_AFTER);
	assert(polar != NULL);
	return (polar);
}

static void
blend_polar(const polar_diag_t *diag, polar_t *tgt_polar,
    const polar_t *src_polar1, double d_alpha)
{
	double blend_tgt_alpha = src_polar1->alpha + d_alpha;
	const polar_t *src_polar2 = find_nearest_polar(diag, blend_tgt_alpha);
	double interp_w, blend_w, Cl, Cd, Cm;
	double alpha_min, alpha_max;

	assert(src_polar2 != NULL);

	alpha_min = MIN(src_polar1->alpha, src_polar1->alpha + d_alpha);
	alpha_max = MAX(src_polar1->alpha, src_polar1->alpha + d_alpha);
	if (tgt_polar->alpha < alpha_min || tgt_polar->alpha > alpha_max)
		return;

	interp_w = iter_fract(tgt_polar->alpha, src_polar1->alpha,
	    src_polar2->alpha, true);
	Cl = wavg(src_polar1->Cl, src_polar2->Cl, interp_w);
	Cd = wavg(src_polar1->Cd, src_polar2->Cd, interp_w);
	Cm = wavg(src_polar1->Cm, src_polar2->Cm, interp_w);

	blend_w = iter_fract(fabs(tgt_polar->alpha - blend_tgt_alpha),
	    0, fabs(d_alpha), true);
	tgt_polar->Cl = wavg(tgt_polar->Cl, Cl, blend_w);
	tgt_polar->Cd = wavg(tgt_polar->Cd, Cd, blend_w);
	tgt_polar->Cm = wavg(tgt_polar->Cm, Cm, blend_w);
}

static void
afl_interp_diag_polar(const polar_diag_t *d1, const polar_diag_t *d2,
    polar_t *p, avl_index_t where)
{
	const polar_t *p1, *p2;
	double w;

	p1 = avl_nearest(&d2->polars, where, AVL_BEFORE);
	p2 = avl_nearest(&d2->polars, where, AVL_AFTER);
	if (p1 == NULL) {
		assert(p2 != NULL);
		blend_polar(d1, p, p2, -blend_range);
	} else if (p2 == NULL) {
		assert(p1 != NULL);
		blend_polar(d1, p, p1, blend_range);
	} else {
		w = iter_fract(p->alpha, p1->alpha, p2->alpha, true);
		p->Cl = wavg(p1->Cl, p2->Cl, w);
		p->Cd = wavg(p1->Cd, p2->Cd, w);
		p->Cm = wavg(p1->Cm, p2->Cm, w);
	}
}

static void
afl_combine_diag(const polar_diag_t *d1, const polar_diag_t *d2)
{
	for (polar_t *p = avl_first(&d1->polars); p != NULL;
	    p = AVL_NEXT(&d1->polars, p)) {
		const polar_t *p2;
		avl_index_t where;

		p2 = avl_find(&d2->polars, p, &where);
		if (p2 != NULL) {
			p->Cl = p2->Cl;
			p->Cd = p2->Cd;
			p->Cm = p2->Cm;
		} else {
			afl_interp_diag_polar(d1, d2, p, where);
		}
	}
}

static void
afl_combine(afl_t *afl, const afl_t *xfoil)
{
	for (polar_diag_t *d1 = avl_first(&afl->diags); d1 != NULL;
	    d1 = AVL_NEXT(&afl->diags, d1)) {
		const polar_diag_t *d2 = avl_find(&xfoil->diags, d1, NULL);

		if (d2 == NULL)
			continue;
		assert(avl_numnodes(&d2->polars) != 0);

		afl_combine_diag(d1, d2);
	}
}

int
main(int argc, const char *argv[])
{
	afl_t *afl;
	afl_t *xfoil;

	if (argc < 3) {
		fprintf(stderr, "Usage: %s <input.afl> <output.afl> "
		    "[polar.txt...]\n", argv[0]);
		return (1);
	}

	afl = afl_parse(argv[1]);
	if (afl == NULL)
		return (1);

	xfoil = afl_alloc();
	for (int i = 3; i < argc; i++) {
		if (!xfoil_polar_parse(xfoil, argv[i]))
			return (1);
	}
	afl_combine(afl, xfoil);
	afl_write(afl, argv[2]);

	afl_free(afl);
	afl_free(xfoil);

	return (0);
}
