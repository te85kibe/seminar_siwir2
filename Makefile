make:
	make -s -f dirichlet.mak clean
	make -s -f dirichlet.mak

neumann:
	make -s -f neumann.mak clean
	make -s -f neumann.mak

seminar:
	make -s -f seminar.mak clean
	make -s -f seminar.mak

clean:
	make -s -f neumann.mak clean
	make -s -f dirichlet.mak clean
