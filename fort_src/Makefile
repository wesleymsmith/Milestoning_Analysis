CC=gfortran

objects = writeindex writerate qelements
all: $(objects)

$(objects): %: %.f
	$(CC) -o $@.exe $<
