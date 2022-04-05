const unsigned long long my_drand48_a(0x5DEECE66D);
const unsigned long long my_drand48_c(0xB);
const unsigned long long my_drand48_m(0x1000000000000);

unsigned long long my_drand48_x(666);

void seed(unsigned long long seed){ my_drand48_x = seed; }

double my_drand48(void)
{
    my_drand48_x = ( my_drand48_a*my_drand48_x + my_drand48_c )%my_drand48_m;
    return (double)my_drand48_x/(double)my_drand48_m;
}
