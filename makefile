run:
	make disease
	./des
	rm des
	make plot
addon:
	make addons
	./add
	rm add
	make plotAddon
	

disease: numerics_disease.o disease.o
	gcc numerics_disease.o disease.o -lgsl -lgslcblas -lm -o des
addons: addon.o numerics_disease.o
	gcc numerics_disease.o addon.o -lgsl -lgslcblas -lm -o add
	
	

disease.o: disease.c
	gcc -c disease.c
addon.o: addon.c
	gcc -Wno-incompatible-pointer-types -c addon.c
	
numerics_disease.o: numerics_disease.c
	gcc -c numerics_disease.c






plot:
	python3 plot.py

plotAddon: 
	python3 plotAddon.py

animate:
	python3 animation.py
	
clean:
	rm *.o
	
