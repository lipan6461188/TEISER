#-O3表示高级优化，-g表示生成警告信息
CFLAGS = -O3 -g

CC = cc

tests : 
	echo $(LFLAGS)

all : mi_find_seed mi_optimize move clean


mi_find_seed : mi_find_seed.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o read_write_motif.o readFASTA.o readicSHAPE.o
	$(CC) $(CFLAGS) -Wall -o mi_find_seed mi_find_seed.c statistics.o dataio.o information.o mi_library.o hashtable.o teiser_functions.o sequences.o matchmaker.o read_write_motif.o readicSHAPE.o readFASTA.o -lm $(LFLAGS)

mi_optimize : mi_optimize.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o readicSHAPE.o
	$(CC) $(CFLAGS) -Wall -o mi_optimize mi_optimize.c teiser_functions.o dataio.o sequences.o hashtable.o matchmaker.o information.o statistics.o read_write_motif.o mi_library.o readFASTA.o readicSHAPE.o -lm $(LFLAGS)

teiser_functions.o : teiser_functions.c teiser_functions.h dataio.o sequences.o hashtable.o matchmaker.o read_write_motif.o information.o statistics.o
	$(CC) $(CFLAGS) -Wall -c teiser_functions.c

readFASTA.o : readFASTA.c readFASTA.h dataio.o sequences.o statistics.o
	$(CC) $(CFLAGS) -Wall -c readFASTA.c

readicSHAPE.o : readicSHAPE.c readicSHAPE.h 
	$(CC) $(CFLAGS) -Wall -c readicSHAPE.c

matchmaker.o : matchmaker.c matchmaker.h
	$(CC) $(CFLAGS) -Wall -c matchmaker.c

read_write_motif.o : read_write_motif.c read_write_motif.h
	$(CC) $(CFLAGS) -Wall -c read_write_motif.c

create_motifs.o : create_motifs.c create_motifs.h
	$(CC) $(CFLAGS) -Wall -c create_motifs.c

mi_library.o : mi_library.c mi_library.h
	$(CC) $(CFLAGS) -Wall -c mi_library.c

statistics.o : statistics.c statistics.h
	$(CC) $(CFLAGS) -Wall -c statistics.c

sequences.o : sequences.c sequences.h
	$(CC) $(CFLAGS) -Wall -c sequences.c

dataio.o : dataio.c dataio.h
	$(CC) $(CFLAGS) -Wall -c dataio.c

information.o : information.c information.h
	$(CC) $(CFLAGS) -Wall -c information.c

hashtable.o : hashtable.c hashtable.h
	$(CC) $(CFLAGS) -Wall -c hashtable.c

move:
	cp ./mi_optimize ./test
	cp ./mi_find_seed ./test

clean:
	rm *.o