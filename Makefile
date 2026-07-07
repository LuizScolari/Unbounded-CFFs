# Unbounded Cover-Free Families — unified build
# Builds a single `cff` binary that dispatches to the RS, Subset and CRS
# constructions. Supports Linux (GCC) and macOS (Clang + Homebrew).

TARGET   = cff
CC       = gcc

SRC_DIR   = src
BUILD_DIR = build

# Operating-system detection
UNAME_S := $(shell uname -s)

# GLib (required by the RS construction)
GLIB_CFLAGS  := $(shell pkg-config --cflags glib-2.0 2>/dev/null)
GLIB_LDFLAGS := $(shell pkg-config --libs   glib-2.0 2>/dev/null)

ifeq ($(UNAME_S),Darwin)
    # macOS: Homebrew + Apple Clang
    BREW_PREFIX := $(shell brew --prefix)
    OMP_CFLAGS   = -Xpreprocessor -fopenmp -I$(BREW_PREFIX)/opt/libomp/include
    OMP_LDFLAGS  = -L$(BREW_PREFIX)/opt/libomp/lib -lomp
    PLATFORM_INCLUDES = -I$(BREW_PREFIX)/include
    PLATFORM_LIBS     = -L$(BREW_PREFIX)/lib
else
    # Linux: GCC
    OMP_CFLAGS   = -fopenmp
    OMP_LDFLAGS  = -fopenmp
    PLATFORM_INCLUDES =
    PLATFORM_LIBS     =
endif

# Compilation / linking flags
COMMON_CFLAGS = -O3 -Wall -Wextra -std=gnu11 -pthread \
                $(PLATFORM_INCLUDES) $(OMP_CFLAGS)

RS_CFLAGS     = $(COMMON_CFLAGS) $(GLIB_CFLAGS) -I$(SRC_DIR)/RS
CRS_CFLAGS    = $(COMMON_CFLAGS) -I$(SRC_DIR)/CRS
SUBSET_CFLAGS = $(COMMON_CFLAGS) -march=native -I$(SRC_DIR)/Subset
MAIN_CFLAGS   = $(COMMON_CFLAGS)

LDFLAGS = $(PLATFORM_LIBS) -lflint -lgmp -lm -pthread \
          $(OMP_LDFLAGS) $(GLIB_LDFLAGS)

# Sources / objects
RS_SRCS     = $(SRC_DIR)/RS/rs_main.c \
              $(SRC_DIR)/RS/cff_builder.c \
              $(SRC_DIR)/RS/cff_file_generator.c
CRS_SRCS    = $(SRC_DIR)/CRS/crs_main.c \
              $(SRC_DIR)/CRS/generate.c \
              $(SRC_DIR)/CRS/Io.c \
              $(SRC_DIR)/CRS/embed.c
SUBSET_SRCS = $(SRC_DIR)/Subset/subset_main.c \
              $(SRC_DIR)/Subset/generate_cff.c \
              $(SRC_DIR)/Subset/cff_io.c
MAIN_SRC    = $(SRC_DIR)/main.c

RS_OBJS     = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(RS_SRCS))
CRS_OBJS    = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(CRS_SRCS))
SUBSET_OBJS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SUBSET_SRCS))
MAIN_OBJ    = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(MAIN_SRC))

OBJS    = $(MAIN_OBJ) $(RS_OBJS) $(CRS_OBJS) $(SUBSET_OBJS)
HEADERS = $(wildcard $(SRC_DIR)/RS/*.h $(SRC_DIR)/CRS/*.h $(SRC_DIR)/Subset/*.h)

# Per-group compilation flags (target-specific variables)
$(RS_OBJS):     CFLAGS := $(RS_CFLAGS)
$(CRS_OBJS):    CFLAGS := $(CRS_CFLAGS)
$(SUBSET_OBJS): CFLAGS := $(SUBSET_CFLAGS)
$(MAIN_OBJ):    CFLAGS := $(MAIN_CFLAGS)

# Rules
.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR) $(TARGET)
