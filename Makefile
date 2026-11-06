CC=cc
CFLAGS=-O2 -flto -Wall
LDFLAGS=-O2 -flto -lm -lopenblas -lgsl -lcxsparse -lklu -lcjson

.PHONY: all clean firmware upload convergence

all: bin/heat_diffusion

convergence: bin/temporal_convergence_test

obj/mesh_1d.o: src/mesh_1d.c
	@mkdir -p obj
	$(CC) -c $(INCLUDE_DIRS) $(CFLAGS) $^ -o $@

obj/fem_1d.o: src/fem_1d.c
	@mkdir -p obj
	$(CC) -c $(INCLUDE_DIRS) $(CFLAGS) $^ -o $@

obj/fem_function_1d.o: src/fem_function_1d.c
	@mkdir -p obj
	$(CC) -c $(INCLUDE_DIRS) $(CFLAGS) $^ -o $@

obj/heat_problem_1d.o: src/heat_problem_1d.c
	@mkdir -p obj
	$(CC) -c $(INCLUDE_DIRS) $(CFLAGS) $^ -o $@

obj/serial_port.o: src/serial_port.c
	@mkdir -p obj
	$(CC) -c $(INCLUDE_DIRS) $(CFLAGS) $^ -o $@

obj/heat_diffusion.o: src/heat_diffusion.c
	@mkdir -p obj
	$(CC) -c $(INCLUDE_DIRS) $(CFLAGS) $^ -o $@

obj/temporal_convergence_test.o: src/temporal_convergence_test.c
	@mkdir -p obj
	$(CC) -c $(INCLUDE_DIRS) $(CFLAGS) $^ -o $@

bin/heat_diffusion: obj/mesh_1d.o obj/fem_1d.o obj/fem_function_1d.o obj/heat_problem_1d.o obj/serial_port.o obj/heat_diffusion.o
	@mkdir -p bin
	$(CC) $(LIB_DIRS) $(LDFLAGS) $^ -o $@

bin/temporal_convergence_test: obj/mesh_1d.o obj/fem_1d.o obj/fem_function_1d.o obj/heat_problem_1d.o obj/temporal_convergence_test.o
	@mkdir -p bin
	$(CC) $(LIB_DIRS) $(LDFLAGS) $^ -o $@

firmware: firmware/firmware.ino
	arduino-cli compile $^

upload: firmware
	arduino-cli upload $^

clean:
	rm -rf obj bin output
