FC = gfortran
FFLAGS = -O2 -Wall
SRC_DIR = src
OBJ_DIR = build

all: extreme

MODULES = utils density bond io linalg
OBJS = \
	$(addprefix $(OBJ_DIR)/,$(addsuffix .o,$(MODULES))) \
	$(OBJ_DIR)/main.o

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	$(FC) $(FFLAGS) -c $< -J $(OBJ_DIR) -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

extreme: $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)

.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) extreme
