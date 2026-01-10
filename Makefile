# ------------------------
# Defaults (override like: make all CASE=qt_rad)
# ------------------------
CASE  ?= no_rad
CASES ?= no_rad qt_rad cld_rad qt_cld_rad

# Tools / entrypoints
JULIA ?= julia
RUN_SH  := ./bin/run_case.sh
POST_SH := ./bin/postprocess.sh
RAD_SCALE ?= 0.001

# ------------------------
# Phony targets
# ------------------------
.PHONY: help run post all run_all post_all all_all clean clean_figures clean_post

help:
	@echo ""
	@echo "Kuang2008 workflow"
	@echo ""
	@echo "Single case:"
	@echo "  make run   CASE=$(CASE)"
	@echo "  make post  CASE=$(CASE)"
	@echo "  make all   CASE=$(CASE)"
	@echo ""
	@echo "All cases:"
	@echo "  make run_all"
	@echo "  make post_all"
	@echo "  make all_all"
	@echo ""
	@echo "Clean:"
	@echo "  make clean_figures CASE=$(CASE)"
	@echo "  make clean_post    CASE=$(CASE)"
	@echo "  make clean         CASE=$(CASE)"
	@echo ""
	@echo "Defaults:"
	@echo "  CASE  = $(CASE)"
	@echo "  CASES = $(CASES)"
	@echo ""

# ------------------------
# Core pipeline
# ------------------------
run:
	$(RUN_SH) $(CASE) $(RAD_SCALE)

post:
	$(POST_SH) $(CASE) $(RAD_SCALE)

all: run post

# ------------------------
# Batch pipeline
# ------------------------
run_all:
	@set -e; \
	for c in $(CASES); do \
		echo "==> run CASE=$$c"; \
		$(RUN_SH) $$c $(RAD_SCALE); \
	done

post_all:
	@set -e; \
	for c in $(CASES); do \
		echo "==> post CASE=$$c"; \
		$(POST_SH) $$c $(RAD_SCALE); \
	done

all_all: run_all post_all

# ------------------------
# Cleaning
# ------------------------
clean_figures:
	rm -rf figures/$(CASE)

clean_post:
	rm -rf output/$(CASE)/post

clean: clean_figures clean_post
