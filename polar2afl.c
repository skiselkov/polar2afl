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

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
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
	struct {
		double	Re;
		double	slope;
		double	intercept;
		double	alpha_min;
		double	alpha_max;
		double	lin_range;
		double	lift_power;
		double	Cl_max;
		double	stall_drop;
		double	stall_power;
		double	stalled_drop;
		double	Cd_min;
		double	min_Cd_Cl;
		double	Cd_alpha10;
		double	Cd_power;
		double	buck_Cl;
		double	buck_width;
		double	buck_depth;
		double	buck_power;
		double	Cm_alpha1;
		double	Cm_alpha2;
		double	Cm_val_alpha_neg20;
		double	Cm_val_alpha1;
		double	Cm_val_alpha2;
		double	Cm_val_alpha_pos20;
	} params;
	avl_tree_t	polars;
	avl_node_t	node;
} polar_diag_t;

typedef struct {
	char		preamble[1024];
	avl_tree_t	diags;
} afl_t;

static bool	do_smooth_polars = true;
static double	edge_blend_range = 10;	/* degrees */

static void afl_free(afl_t *afl);
static polar_diag_t *smooth_diag(polar_diag_t *d);

ssize_t
my_getline(char **line_p, size_t *cap_p, FILE *fp)
{
	char *line;
	size_t cap, n = 0;

	assert(line_p != NULL);
	line = *line_p;
	assert(cap_p != NULL);
	cap = *cap_p;
	assert(fp != NULL);

	do {
		if (n + 1 >= cap) {
			cap += 256;
			line = realloc(line, cap);
		}
		assert(n < cap);
		if (fgets(&line[n], cap - n, fp) == NULL) {
			if (n != 0)
				break;
			*line_p = line;
			*cap_p = cap;
			return (-1);
		}
		n = strlen(line);
	} while (n > 0 && line[n - 1] != '\n');

	*line_p = line;
	*cap_p = cap;

	return (n);
}

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

	if (pa->params.Re < pb->params.Re)
		return (-1);
	if (pa->params.Re > pb->params.Re)
		return (1);
	return (0);
}

static int
polar_compar(const void *a, const void *b)
{
	const polar_t *pa = a, *pb = b;
	double alpha_a = round(pa->alpha * 10.0) / 10.0;
	double alpha_b = round(pb->alpha * 10.0) / 10.0;

	if (alpha_a < alpha_b)
		return (-1);
	if (alpha_a > alpha_b)
		return (1);
	return (0);
}

static void
polar_diag_free(polar_diag_t *diag)
{
	void *cookie = NULL;
	polar_t *polar;

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

	n = my_getline(&line, &linecap, fp);
	/* Probably an EOF, just exit */
	if (n <= 0)
		goto errout;
	if (sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
	    "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
	    &diag->params.Re, &diag->params.slope, &diag->params.intercept,
	    &diag->params.alpha_min, &diag->params.alpha_max,
	    &diag->params.lin_range, &diag->params.lift_power,
	    &diag->params.Cl_max, &diag->params.stall_drop,
	    &diag->params.stall_power, &diag->params.stalled_drop,
	    &diag->params.Cd_min, &diag->params.min_Cd_Cl,
	    &diag->params.Cd_alpha10, &diag->params.Cd_power,
	    &diag->params.buck_Cl, &diag->params.buck_width,
	    &diag->params.buck_depth, &diag->params.buck_power,
	    &diag->params.Cm_alpha1, &diag->params.Cm_alpha2,
	    &diag->params.Cm_val_alpha_neg20,
	    &diag->params.Cm_val_alpha1, &diag->params.Cm_val_alpha2,
	    &diag->params.Cm_val_alpha_pos20) != 25) {
		fprintf(stderr, "%s:%d: error parsing polar parameters\n",
		    filename, *line_nr);
		goto errout;
	}
	diag->params.Re *= 1000000.0;
	(*line_nr)++;
	/* skip the "alpha cl cd cm:" line */
	if (my_getline(&line, &linecap, fp) <= 0) {
		fprintf(stderr, "%s:%d: premature EOF\n", filename, *line_nr);
		goto errout;
	}
	(*line_nr)++;

	for (; my_getline(&line, &linecap, fp) > 0; (*line_nr)++) {
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

	for (line_nr = 0; my_getline(&line, &linecap, fp) > 0; line_nr++) {
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
			    "Re = %f\n", filename, diag->params.Re);
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

	fprintf(fp, "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f "
	    "%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f "
	    "%.5f %.5f %.5f\n",
	    diag->params.Re / 1000000, diag->params.slope,
	    diag->params.intercept, diag->params.alpha_min,
	    diag->params.alpha_max, diag->params.lin_range,
	    diag->params.lift_power, diag->params.Cl_max,
	    diag->params.stall_drop, diag->params.stall_power,
	    diag->params.stalled_drop, diag->params.Cd_min,
	    diag->params.min_Cd_Cl, diag->params.Cd_alpha10,
	    diag->params.Cd_power, diag->params.buck_Cl,
	    diag->params.buck_width, diag->params.buck_depth,
	    diag->params.buck_power, diag->params.Cm_alpha1,
	    diag->params.Cm_alpha2, diag->params.Cm_val_alpha_neg20,
	    diag->params.Cm_val_alpha1, diag->params.Cm_val_alpha2,
	    diag->params.Cm_val_alpha_pos20);
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

static double
lerp(double x, double x1, double y1, double x2, double y2)
{
	assert(x1 != x2);
	return (((x - x1) / (x2 - x1)) * (y2 - y1) + y1);
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

	while (my_getline(&line, &linecap, fp) > 0) {
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
				diag->params.Re = mantissa * pow(10, exponent);
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
	/*
	 * Pass over the entire diagram and do a linear interpolation
	 * of any missing pieces.
	 */
	for (polar_t *pt = avl_first(&diag->polars), *pt_next = NULL;
	    pt != NULL; pt = pt_next) {
		pt_next = AVL_NEXT(&diag->polars, pt);
		if (pt_next == NULL)
			break;
		for (double alpha = pt->alpha + 0.1;
		    alpha + 0.01 < pt_next->alpha;
		    alpha = (round(alpha * 10.0) + 1.0) / 10.0) {
			polar_t *polar = safe_calloc(1, sizeof (*polar));

			polar->alpha = alpha;
			polar->Cl = lerp(alpha, pt->alpha, pt->Cl,
			    pt_next->alpha, pt_next->Cl);
			polar->Cd = lerp(alpha, pt->alpha, pt->Cd,
			    pt_next->alpha, pt_next->Cd);
			polar->Cm = lerp(alpha, pt->alpha, pt->Cm,
			    pt_next->alpha, pt_next->Cm);
			avl_add(&diag->polars, polar);
		}
	}

	if (do_smooth_polars)
		diag = smooth_diag(diag);
	if (avl_find(&xfoil->diags, diag, &where) != NULL) {
		fprintf(stderr, "%s: is a duplicate of an existing polar "
		    "at Re = %fM\n", filename, diag->params.Re / 1000000);
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
		blend_polar(d1, p, p2, -edge_blend_range);
	} else if (p2 == NULL) {
		assert(p1 != NULL);
		blend_polar(d1, p, p1, edge_blend_range);
	} else {
		w = iter_fract(p->alpha, p1->alpha, p2->alpha, true);
		p->Cl = wavg(p1->Cl, p2->Cl, w);
		p->Cd = wavg(p1->Cd, p2->Cd, w);
		p->Cm = wavg(p1->Cm, p2->Cm, w);
	}
}

/*
 * Uses a Gaussian interpolator to smooth out a polar diagram. The passed
 * diagram is then freed and the smoothed-out version is returned. The
 * should substitute the original diagram with the newly smoothed-out one.
 */
static polar_diag_t *
smooth_diag(polar_diag_t *d1)
{
	polar_diag_t *d2;
	const double kernel[] = { 0.06136, 0.24477, 0.38774, 0.24477, 0.06136 };

	d2 = polar_diag_alloc();
	memcpy(&d2->params, &d1->params, sizeof (d1->params));

	for (polar_t *p = avl_first(&d1->polars); p != NULL;
	    p = AVL_NEXT(&d1->polars, p)) {
		polar_t *p_l, *p_ll, *p_r, *p_rr;
		polar_t *np = safe_calloc(1, sizeof (*np));

		memcpy(np, p, sizeof (*np));
		avl_add(&d2->polars, np);

		p_l = AVL_PREV(&d1->polars, p);
		p_r = AVL_NEXT(&d1->polars, p);
		if (p_l == NULL || p_r == NULL)
			continue;
		p_ll = AVL_PREV(&d1->polars, p_l);
		p_rr = AVL_NEXT(&d1->polars, p_r);
		if (p_ll == NULL || p_rr == NULL)
			continue;

		np->Cl = (p_ll->Cl * kernel[0] + p_l->Cl * kernel[1] +
		    p->Cl * kernel[2] + p_r->Cl * kernel[3] +
		    p_rr->Cl * kernel[4]);
		np->Cd = (p_ll->Cd * kernel[0] + p_l->Cd * kernel[1] +
		    p->Cd * kernel[2] + p_r->Cd * kernel[3] +
		    p_rr->Cd * kernel[4]);
		np->Cm = (p_ll->Cm * kernel[0] + p_l->Cm * kernel[1] +
		    p->Cm * kernel[2] + p_r->Cm * kernel[3] +
		    p_rr->Cm * kernel[4]);
	}

	polar_diag_free(d1);

	return (d2);
}

static void
afl_combine_diag(polar_diag_t *d1, const polar_diag_t *d2)
{
	double Cl_max = -INFINITY;
	double Cl_max_alpha = NAN;
	double Cl_min = INFINITY;
	double Cl_min_alpha = NAN;
	double Cd_min = INFINITY;
	double Cd_min_Cl = NAN;

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

		if (p->Cl > Cl_max) {
			Cl_max = p->Cl;
			Cl_max_alpha = p->alpha;
		}
		if (p->Cl < Cl_min) {
			Cl_min = p->Cl;
			Cl_min_alpha = p->alpha;
		}
		if (p->Cd < Cd_min) {
			Cd_min = p->Cd;
			Cd_min_Cl = p->Cl;
		}
	}

	assert(isfinite(Cl_max));
	assert(!isnan(Cl_max_alpha));
	assert(isfinite(Cl_min));
	assert(!isnan(Cl_min_alpha));
	assert(isfinite(Cd_min));
	assert(!isnan(Cd_min_Cl));

	d1->params.Cl_max = Cl_max;
	d1->params.alpha_max = clamp(Cl_max_alpha, -19.999, 19.999);
	d1->params.alpha_min = clamp(Cl_min_alpha, -19.999, 19.999);
	d1->params.Cd_min = Cd_min;
	d1->params.min_Cd_Cl = Cd_min_Cl;
}

static bool
afl_combine(afl_t *afl, const afl_t *xfoil)
{
	bool result = true;

	for (polar_diag_t *d_xfoil = avl_first(&xfoil->diags); d_xfoil != NULL;
	    d_xfoil = AVL_NEXT(&afl->diags, d_xfoil)) {
		const polar_diag_t *d_afl =
		    avl_find(&afl->diags, d_xfoil, NULL);

		if (d_afl == NULL) {
			fprintf(stderr, "ERROR: XFoil polar Re=%.3fM has no "
			    "matching polar in the AFL file.\n"
			    "  Open the input AFL file in Airfoil Maker and "
			    "add the missing Re number first,\n"
			    "  then re-run polar2afl.\n",
			    d_xfoil->params.Re / 1000000.0);
			result = false;
		}
	}

	for (polar_diag_t *d_afl = avl_first(&afl->diags); d_afl != NULL;
	    d_afl = AVL_NEXT(&afl->diags, d_afl)) {
		const polar_diag_t *d_xfoil =
		    avl_find(&xfoil->diags, d_afl, NULL);

		if (d_xfoil == NULL) {
			fprintf(stderr, "ERROR: polar Re=%.3fM in AFL file "
			    "has no matching XFoil polar,\n"
			    "  leaving unmodified.\n",
			    d_afl->params.Re / 1000000.0);
			result = false;
			continue;
		}
		assert(avl_numnodes(&d_xfoil->polars) != 0);

		afl_combine_diag(d_afl, d_xfoil);
	}

	return (result);
}

static void
print_usage(FILE *fp, const char *progname)
{
	fprintf(fp, "Usage: %s [-hs] [-e <range>] <input.afl> <output.afl> "
	    "[polar.txt...]\n"
	    "\n"
	    "  Copyright 2019 Saso Kiselkov. All rights reserved.\n"
	    "\n"
	    "  This utility takes an input .afl file produced by X-Plane's\n"
	    "  Airfoil Maker and a number of polar graphs from XFoil and\n"
	    "  splices the XFoil data into the .afl file, generating a new\n"
	    "  output .afl file in the process. If you are using  XFLR5 to\n"
	    "  generate the polars, you can batch export all polars by\n"
	    "  selecting Polars -> Export all -> to text format.\n"
	    "  Please note: before attempt to splice the XFoil polars an .afl\n"
	    "  file, make sure that the file already contains all the Re\n"
	    "  numbers you want to insert. For example, if you've generated\n"
	    "  XFoil polars for Re=3.5M, 4.5M and 5.5M, make sure the\n"
	    "  corresponding Re tabs are already present in the .afl file.\n"
	    "  Then save the file in Airfoil Maker and use polar2afl.\n"
	    "  polar2afl will NOT splice in any polars for Re numbers that\n"
	    "  are not already present in the .afl file.\n"
	    "\n"
	    "Example:\n"
	    "  Take two polars named `NACA2412_Re2.5.txt' and "
	    "`NACA2412_Re3.5.txt'\n"
	    "  in the current directory and splice them into an airfoil file\n"
	    "  named `foo.afl', writing the result to a new file named "
	    "`bar.afl'.\n"
	    "\n"
	    "  $ %s foo.afl bar.afl NACA2412_Re2.5.txt NACA2412_Re3.5.txt\n"
	    "\n"
	    "Options:\n"
	    " -h: Shows this help screen and exits.\n"
	    " -s: DON'T smooth polar diagrams. The polars produced by XFoil\n"
	    "     sometimes contain jagged jumps, especially when the solver\n"
	    "     couldn't converge on some points. To avoid such weirdness\n"
	    "     propagating into the .afl file, we do a guassian smoothing\n"
	    "     pass to get rid of such artifacts. Passing the '-s' option\n"
	    "     DISABLES this smoothing pass.\n"
	    " -e <range>: selects the alpha range over which we interpolate\n"
	    "     between the XFoil polars and the original input .afl file.\n"
	    "     The edges of the XFoil polars generally don't smoothly\n"
	    "     attach to the original .afl file, so over a range of alphas\n"
	    "     outside of the polar data range, we gradually interpolate\n"
	    "     between the original .afl curve and our polar data curve.\n"
	    "     The default is to interpolate the edge fit over 10 degrees\n"
	    "     of alpha range.\n",
	    progname, progname);
}

int
main(int argc, char *argv[])
{
	afl_t *afl;
	afl_t *xfoil;
	int opt;
	bool result;

	while ((opt = getopt(argc, argv, "hse:")) != -1) {
		switch (opt) {
		case 'h':
			print_usage(stdout, argv[0]);
			return (0);
		case 's':
			do_smooth_polars = false;
			break;
		case 'e':
			edge_blend_range = atof(optarg);
			if (edge_blend_range < 0.1) {
				fprintf(stderr, "Invalid edge blend range, "
				    "must be greater than or equal to 0.1\n");
				return (1);
			}
			break;
		default:
			print_usage(stderr, argv[0]);
			return (1);
		}
	}

	argc -= optind - 1;
	argv += optind - 1;

	if (argc < 3) {
		fprintf(stderr, "Missing arguments. Try \"%s -h\" for help.\n",
		    argv[0]);
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
	result = afl_combine(afl, xfoil);
	afl_write(afl, argv[2]);

	afl_free(afl);
	afl_free(xfoil);

	return (result ? 0 : 1);
}
