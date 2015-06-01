make:
	make -s -f seminar.mak clean
	make -s -f seminar.mak

clean:
	make -s -f neumann.mak clean
	make -s -f dirichlet.mak clean
