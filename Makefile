all:
	cd src; make

clean: 
	cd src; make clean

cleanall: 
	rm -rf output/
	rm -rf obj/
