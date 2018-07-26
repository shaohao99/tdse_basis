include ../common.mk

all: libbasis.a

LIBSRCDIR=../$(LIBSRCDIRNAME)
LIBBASISDIR=../$(LIBBASIS)

libbasis.a: $(LIBBASISDIR)/libbasis.a
	ln $^ $@

#$(LIBSRCDIR)/libnpsf.a:
#	$(MAKE) -C $(LIBSRCDIR)

clean: 
	@rm -f *.a

.PHONY: all clean
