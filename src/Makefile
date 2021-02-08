CXX=g++
CXXFLAGS=  -g -O2
BIN := ../bin
LDFLAGS=-lz
INCLUDE=  -L/usr/lib64/   -L/usr/lib/ -L./include/zlib/
all: $(BIN)/LDBlockShow

$(BIN)/LDBlockShow: $(BIN)/LDBlockShow.o 
	$(CXX)   $^ -o $@   $(LDFLAGS)  $(INCLUDE) 

$(BIN)/%.o: %.cpp
	$(CXX)  -c $(CXXFLAGS) $< -o $@    $(INCLUDE)
	chmod  755 ${BIN}/*

clean:
	$(RM) $(BIN)/*.o
