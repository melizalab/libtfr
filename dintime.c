
/*  copyright 1990  Richard H. Krukar all rights reserved

	Permission granted to buy, sell, or steal this software is granted.
    The author retains the right to distribute this software freely, but
    is not responsible for it's quality or maintainance. */

/*  A split radix ( 2/4 ) FFT based on ( reference here )  This routine
is meant to be fairly efficient, but not exeptionally.  Recursion is
used to sort out the butterflys.  This should not matter because as
soon as the length reaches UNROL, everything is inlined elsewhere */

struct complex { float r,i; };

#include "fft.h"

extern void dint(float *x, int length, float *wtab);

void dintime(struct complex *x, int length, float *wtab )
{
    int qlen, step1, step2, c1, c2, s1, s2;
    struct complex *ptr1, *ptr2, *ptr3, *ptr4, tmp;

    if(length<= UNROL ) {
	dint( (float *)x, length, wtab);
	return; }

    s1 = s2 = TLEN/4;
    c1 = c2 = 0;
    step1 = TLEN/length;
    step2 = 3*step1;

    qlen = length/4;
    for(ptr1 = x; ptr1 < &x[qlen]; ptr1++) {
	ptr2 = ptr1+qlen;
	ptr3 = ptr2+qlen;
	ptr4 = ptr3+qlen;

	tmp.r = ptr2->r - ptr4->r; /* First Butterfly */
	tmp.i = ptr2->i - ptr4->i;
	ptr2->r += ptr4->r;
	ptr2->i += ptr4->i;

	ptr4->r = tmp.i; ptr4->i = -tmp.r; /* mult by -j */

	tmp.r = ptr1->r - ptr3->r; /* Second Butterfly */
	tmp.i = ptr1->i - ptr3->i;
	ptr1->r += ptr3->r;
	ptr1->i += ptr3->i;

	ptr3->r = tmp.r + ptr4->r; /* Third Butterfly */
	ptr3->i = tmp.i + ptr4->i;
	ptr4->r = tmp.r - ptr4->r;
	ptr4->i = tmp.i - ptr4->i;
	
/* Now multiply ptr3 and ptr4 by the appropriate weights */

	tmp.r = ptr3->r*wtab[c1] - ptr3->i*wtab[s1];
	ptr3->i = ptr3->r*wtab[s1] + ptr3->i*wtab[c1];
	ptr3->r = tmp.r;

	tmp.r = ptr4->r*wtab[c2] - ptr4->i*wtab[s2];
	ptr4->i = ptr4->r*wtab[s2] + ptr4->i*wtab[c2];
	ptr4->r = tmp.r;

	c1 += step1; s1 += step1;
	c2 += step2; s2 += step2;
    }
    dintime(&x[0], (length/2), wtab);
    dintime(&x[(length/2)], (length/4), wtab);
    dintime(&x[((length*3)/4)], (length/4), wtab);
}

