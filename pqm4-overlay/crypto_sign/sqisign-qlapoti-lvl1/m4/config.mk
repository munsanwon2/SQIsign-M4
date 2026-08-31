# pqm4 target-specific options for crypto_sign/sqisign-qlapoti-lvl1/m4.
# mupq/mk/schemes.mk loads this file automatically.

SQIQL_M4_IMPL := crypto_sign_sqisign-qlapoti-lvl1_m4
SQIQL_M4_DIR := $(patsubst %/,%,$(dir $(lastword $(MAKEFILE_LIST))))
SQIQL_M4_BACKEND_DIRS := $(shell if test -d "$(SQIQL_M4_DIR)/backend"; then find "$(SQIQL_M4_DIR)/backend" -type d -print; fi)
SQIQL_M4_INCLUDE_FLAGS := -I$(SQIQL_M4_DIR) $(addprefix -I,$(SQIQL_M4_BACKEND_DIRS))
SQIQL_M4_AIO_FLAGS := -Os -Wvla -Wframe-larger-than=2048 -fstack-usage -std=gnu11
SQIQL_M4_BACKEND_AUDIT_FLAGS := $(SQIQL_M4_AIO_FLAGS) -Werror=vla

ifeq ($(PLATFORM),stm32f4discovery)
elf/crypto_sign_sqisign-qlapoti-lvl1_m4_%.elf: LDSCRIPT = ldscripts/stm32f4discovery_sqisign_qlapoti.ld
SQIQL_M4_SECTION_FLAG := -DQLW_SECTION_NAME=\".ccmram.qlapoti_workspace\" -DSQISIGN_M4_F407_BANKED_RAM=1
else
SQIQL_M4_SECTION_FLAG := -DQLW_SECTION_NAME=\".bss.qlapoti_workspace\"
endif

# Build the implementation archive as a backend-only compile audit.  pqm4 AIO
# targets compile their harness and the implementation in one compiler command,
# and testvectors.c intentionally contains a VLA.  Keeping -Werror=vla on the
# archive preserves the backend gate without rejecting the upstream harness.
obj/lib$(SQIQL_M4_IMPL).a: CPPFLAGS += $(SQIQL_M4_INCLUDE_FLAGS) $(SQIQL_M4_SECTION_FLAG)
obj/lib$(SQIQL_M4_IMPL).a: CFLAGS += $(SQIQL_M4_INCLUDE_FLAGS) $(SQIQL_M4_SECTION_FLAG) $(SQIQL_M4_BACKEND_AUDIT_FLAGS)
elf/$(SQIQL_M4_IMPL)_%.elf: CPPFLAGS += $(SQIQL_M4_INCLUDE_FLAGS) $(SQIQL_M4_SECTION_FLAG)
elf/$(SQIQL_M4_IMPL)_%.elf: CFLAGS += $(SQIQL_M4_INCLUDE_FLAGS) $(SQIQL_M4_SECTION_FLAG) $(SQIQL_M4_AIO_FLAGS)
# The curated backend is its own audited FIPS202 provider.  pqm4 AIO normally
# compiles mupq/common/fips202.c into every firmware; retain all other common
# sources while excluding that one duplicate provider for these ELFs.  Filter
# SYMCRYPTO_SRC (not LIBDEPS) so pqm4's later per-harness NO_RANDOMBYTES
# expansion remains live, especially for deterministic testvectors.
elf/$(SQIQL_M4_IMPL)_%.elf: SYMCRYPTO_SRC := $(filter-out mupq/common/fips202.c,$(SYMCRYPTO_SRC))
elf/$(SQIQL_M4_IMPL)_%.elf: LINKDEPS += obj/lib$(SQIQL_M4_IMPL).a
