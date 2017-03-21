
FC = gfortran -cpp

fno:
	$(FC) -o fno funspeed.f90

f-O:
	$(FC) -o f-O -O funspeed.f90

f-O3:
	$(FC) -o f-O3 -O3 funspeed.f90

ffast:
	$(FC) -o ffast -Ofast -O3 funspeed.f90

fd-no:
	$(FC) -DDOUBLE_CLOCK -o fd-no funspeed.f90

fd-O:
	$(FC) -DDOUBLE_CLOCK -o fd-O -O funspeed.f90

fd-O3:
	$(FC) -DDOUBLE_CLOCK  -o fd-O3 -O3 funspeed.f90

fd-fast:
	$(FC) -DDOUBLE_CLOCK -o fd-fast -Ofast -O3 funspeed.f90


frt-O:
	$(FC) -DDOUBLE_CLOCK -lrt -o frt-O -O funspeed.f90

frt-O3:
	$(FC) -DDOUBLE_CLOCK -lrt -o frt-O3 -O3 funspeed.f90


frt-fast:
	$(FC) -DDOUBLE_CLOCK -lrt -o frt-fast -Ofast -O3 funspeed.f90


run: f-O fno f-O3 ffast fd-no fd-O fd-O3 fd-fast frt-O frt-O3 frt-fast
	echo "No optimisation"; ./fno ;\
	echo "Min optimisation -O"; ./f-O ;\
	echo "Level 3 optimisation '-O3'"; ./f-O3 ;\
	echo "Fast optimisation '-Ofast'"; ./ffast ;\
	echo "DC No optimisation"; ./fd-no ;\
	echo "DC Min optimisation -O"; ./fd-O ;\
	echo "DC Level 3 optimisation '-O3'"; ./fd-O3 ;\
	echo "DC Fast optimisation '-Ofast'"; ./fd-fast ;\
	echo "RT Min optimisation "; ./frt-O ;\
	echo "RT Level 3 optimisation '-O3 -lrt'"; ./frt-O3 ;\
	echo "RT Fast otimisation '-Ofast -lrt'"; ./frt-fast ;\


clean:
	rm -f f-O fno f-O3 ffast fd-no fd-O fd-O3 fd-fast frt-O frt-O3 frt-fast

all: clean run
