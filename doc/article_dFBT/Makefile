FNAME=FJTrafo_v3
#FNAME=diff
#FNAME=FJTrafo_v1

ps: dvi
	dvips $(FNAME).dvi
	
dvi:
	latex $(FNAME)
	bibtex $(FNAME)
	latex $(FNAME)
	latex $(FNAME)


pdf: ps
	ps2pdf $(FNAME).ps

view: ps
	open $(FNAME).ps &

clean:
	rm -f $(FNAME).{spl,bbl,ps,log,aux,out,dvi,blg}
	rm -f *~

.PHONY: all clean
