
/*  copyright 1990  Richard H. Krukar all rights reserved

	Permission granted to buy, sell, or steal this software is granted.
    The author retains the right to distribute this software freely, but
    is not responsible for it's quality or maintainance. */

void bitrev512( float x[], int length);

void bitrev( float *x, int length)
{
	int i, irev, j, bits;
	float temp, *ptr1, *ptr2;

	if(length <= 512) {
	    bitrev512(x, length);
	    return; }

	for(i=0;i<length;i++) {
	    irev = 0;
	    for(j=1, bits=(length>>1);j<length;(j=j<<1),(bits=bits>>1))
		if(i&j) irev |= bits;
	    if(i < irev) {
		ptr1 = &(x[i<<1]); ptr2 = &(x[irev<<1]); 
		temp= *ptr1; *ptr1 = *ptr2; *ptr2 = temp;
		ptr1++; ptr2++;
		temp= *ptr1; *ptr1 = *ptr2; *ptr2 = temp;
	    }
	}
}
