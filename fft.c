
/*  copyright 1990  Richard H. Krukar all rights reserved

	Permission granted to buy, sell, or steal this software is granted.
    The author retains the right to distribute this software freely, but
    is not responsible for it's quality or maintainance. */

#include "fft.h"
struct complex { float r,i; };

extern float wtab[];
extern void dintime(struct complex *x, int length, float *wtab );
extern void idintime( struct complex x[], int length, float *wtab );
extern void bitrev( float *x, int length);
extern void bitrev512( float x[], int length);

void fft(float x[], int length )
{
    dintime( (struct complex *)x, length, wtab );
    bitrev( x, length );
}


void ifft( float x[], int length )
{
    float *fptr, *eptr, scale;

    idintime( (struct complex *)x, length, wtab );
    bitrev( x, length );

    scale = (1.0/length); eptr = &x[length<<1];
    for(fptr = &x[0]; fptr<eptr; fptr++)
	*fptr *= scale;
}

void rfft( struct complex x[], int length )
{
    int hlen, step, i, c, s;
    struct complex odd, even;

    hlen = length/2;

    dintime( x, hlen, wtab );
    bitrev( (float *)x, hlen );

    c = step = TLEN/length;
    s = c + TLEN/4;

    for(i=1; i<=(length/4); i++) {
	even.r = x[i].r + x[hlen-i].r;
	even.i = x[i].i - x[hlen-i].i;

	odd.r = x[i].r - x[hlen-i].r;
	odd.i = x[i].i + x[hlen-i].i;

	x[i].i = odd.i * wtab[s] - odd.r * wtab[c]; /* use x as a tmp */
	odd.r = odd.r * wtab[s] + odd.i * wtab[c];
	odd.i = x[i].i;

	x[i].r = .5*(even.r+odd.r);
	x[i].i = .5*(even.i+odd.i);
	x[hlen-i].r = .5*(even.r-odd.r);
	x[hlen-i].i = .5*(odd.i-even.i);
	c += step;
	s += step;
	if( s > TLEN ) s -= TLEN;
    }
    even.r = x[0].r + x[0].i;
    x[0].i = x[0].r - x[0].i;
    x[0].r = even.r;
}


void rifft( struct complex x[], int length )
{
    int hlen, step, i, c, s;
    struct complex odd, even;
    float *fptr, *eptr, scale;

    hlen = length/2;
    c = step = TLEN/length;
    s = c + 3*TLEN/4;

    for(i=1; i<=(length/4); i++) {
	even.r = x[i].r + x[hlen-i].r;
	even.i = x[i].i - x[hlen-i].i;

	odd.r = x[i].r - x[hlen-i].r;
	odd.i = x[i].i + x[hlen-i].i;

	x[i].i = odd.i * wtab[s] - odd.r * wtab[c]; /* use x as a tmp */
	odd.r = odd.r * wtab[s] + odd.i * wtab[c];
	odd.i = x[i].i;

	x[i].r = .5*(even.r-odd.r);
	x[i].i = .5*(even.i-odd.i);
	x[hlen-i].r = .5*(even.r+odd.r);
	x[hlen-i].i = .5*(-odd.i-even.i);
	c += step;
	s += step;
	if( s > TLEN ) s -= TLEN;
    }
    even.r = .5*(x[0].r + x[0].i);
    x[0].i = .5*(x[0].r - x[0].i);
    x[0].r = even.r;
    idintime( x, hlen, wtab );
    bitrev( (float *)x, hlen );

    scale = (1.0/hlen); eptr = (float *)&x[length>>1];
    for(fptr = (float *)&x[0]; fptr<eptr; fptr++)
	*fptr *= scale;
}
