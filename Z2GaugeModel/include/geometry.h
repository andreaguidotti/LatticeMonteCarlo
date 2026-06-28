
#ifndef GEOMETRY_H
#define GEOMETRY_H

inline long int dirgeo(long int lex, int i, long int volume)
{
    return lex + volume * i;
}

void cart_to_lex(long int *lex, int *cartcoord, int dim, int L);

void lex_to_cart(long int lex,  int *cartcoord, int dim, int L);

void init_geo(long int *  restrict nnp, long int * restrict nnm, int L, int dim);

#endif
