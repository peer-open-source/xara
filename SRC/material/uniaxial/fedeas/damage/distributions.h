/* Provided under the Zero-Clause BSD (0BSD) License 
 *
 * Copyright (C) 2021 Claudio Perez
 *
 * Permission to use, copy, modify, and/or distribute this software 
 * for any purpose with or without fee is hereby granted.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL
 * THE AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR
 * CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM
 * LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,
 * NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN
 * CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */


#ifdef __cplusplus
extern "C" {
#endif

typedef double (cdffun)(double, double, double);
extern cdffun logncdf;
extern cdffun lognpdf;

/* Weibull PDF
 *
 * f(x|a,b) = (b/a)(x/a)^(b-1) exp(-(x/a)^b)
 */
double wblpdf(double x, double a, double b);
double wblcdf(double x, double a, double b);


/* Lognormal PDF
 */
double lognpdf(double x, double mean, double stddev);
double logncdf(double x, double mean, double stddev);


/* PDF for the 2-parameter Beta distribution.
 */
double betapdf(double x, double a, double b);

/* CDF for the 2-parameter Beta distribution.
 *
 * F(x; a,b) = B(x; a,b) / B(a, b) = I_x(a,b) */
double betacdf(double x, double a, double b);

#ifdef __cplusplus
}
#endif