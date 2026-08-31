# pqm4 options for the portable host/reference companion under
# mupq/crypto_sign/sqisign-qlapoti-lvl1/ref.

SQIQL_REF_IMPL := mupq_crypto_sign_sqisign-qlapoti-lvl1_ref
SQIQL_REF_DIR := $(patsubst %/,%,$(dir $(lastword $(MAKEFILE_LIST))))
SQIQL_REF_BACKEND_DIRS := $(shell if test -d "$(SQIQL_REF_DIR)/backend"; then find "$(SQIQL_REF_DIR)/backend" -type d -print; fi)
SQIQL_REF_INCLUDE_FLAGS := -I$(SQIQL_REF_DIR) $(addprefix -I,$(SQIQL_REF_BACKEND_DIRS))
SQIQL_REF_AIO_FLAGS := -Os -Wvla -Wframe-larger-than=2048 -fstack-usage -std=gnu11
SQIQL_REF_BACKEND_AUDIT_FLAGS := $(SQIQL_REF_AIO_FLAGS) -Werror=vla
SQIQL_REF_HOST_AIO_FLAGS := -Wvla -Wframe-larger-than=2048 -std=gnu11
SQIQL_REF_HOST_BACKEND_AUDIT_FLAGS := $(SQIQL_REF_HOST_AIO_FLAGS) -Werror=vla
SQIQL_REF_HOST_AUDIT_LIB := obj-host/lib$(SQIQL_REF_IMPL)_vla_audit.a

ifeq ($(PLATFORM),stm32f4discovery)
elf/mupq_crypto_sign_sqisign-qlapoti-lvl1_ref_%.elf: LDSCRIPT = ldscripts/stm32f4discovery_sqisign_qlapoti.ld
SQIQL_REF_TARGET_SECTION_FLAG := -DQLW_SECTION_NAME=\".ccmram.qlapoti_workspace\" -DSQISIGN_M4_F407_BANKED_RAM=1
else
SQIQL_REF_TARGET_SECTION_FLAG := -DQLW_SECTION_NAME=\".bss.qlapoti_workspace\"
endif
# The curated library excludes its hosted OS RNG.  Opt into only its named
# platform-entropy hook; do not enable TARGET_ARM ABI/layout behavior in this
# portable host build.  The harness supplies deterministic PQCLEAN_randombytes.
SQIQL_REF_HOST_SECTION_FLAG := -DQLW_SECTION_NAME=\".bss.qlapoti_workspace\" -DSQISIGN_USE_PLATFORM_RANDOMBYTES=1

obj/lib$(SQIQL_REF_IMPL).a: CPPFLAGS += $(SQIQL_REF_INCLUDE_FLAGS) $(SQIQL_REF_TARGET_SECTION_FLAG)
obj/lib$(SQIQL_REF_IMPL).a: CFLAGS += $(SQIQL_REF_INCLUDE_FLAGS) $(SQIQL_REF_TARGET_SECTION_FLAG) $(SQIQL_REF_BACKEND_AUDIT_FLAGS)
elf/$(SQIQL_REF_IMPL)_%.elf: CPPFLAGS += $(SQIQL_REF_INCLUDE_FLAGS) $(SQIQL_REF_TARGET_SECTION_FLAG)
elf/$(SQIQL_REF_IMPL)_%.elf: CFLAGS += $(SQIQL_REF_INCLUDE_FLAGS) $(SQIQL_REF_TARGET_SECTION_FLAG) $(SQIQL_REF_AIO_FLAGS)
elf/$(SQIQL_REF_IMPL)_%.elf: SYMCRYPTO_SRC := $(filter-out mupq/common/fips202.c,$(SYMCRYPTO_SRC))
elf/$(SQIQL_REF_IMPL)_%.elf: LINKDEPS += obj/lib$(SQIQL_REF_IMPL).a

# Compile the implementation alone with -Werror=vla before linking pqm4's host
# testvectors harness, whose print buffer intentionally is a VLA.
$(SQIQL_REF_HOST_AUDIT_LIB): $(call hostobjs,$(call schemesrc,$(SQIQL_REF_DIR)))
$(SQIQL_REF_HOST_AUDIT_LIB): HOST_CPPFLAGS += $(SQIQL_REF_INCLUDE_FLAGS) $(SQIQL_REF_HOST_SECTION_FLAG)
$(SQIQL_REF_HOST_AUDIT_LIB): HOST_CFLAGS += $(SQIQL_REF_INCLUDE_FLAGS) $(SQIQL_REF_HOST_SECTION_FLAG) $(SQIQL_REF_HOST_BACKEND_AUDIT_FLAGS)
bin-host/$(SQIQL_REF_IMPL)_testvectors: HOST_CPPFLAGS += $(SQIQL_REF_INCLUDE_FLAGS) $(SQIQL_REF_HOST_SECTION_FLAG)
bin-host/$(SQIQL_REF_IMPL)_testvectors: HOST_CFLAGS += $(SQIQL_REF_INCLUDE_FLAGS) $(SQIQL_REF_HOST_SECTION_FLAG) $(SQIQL_REF_HOST_AIO_FLAGS)
bin-host/$(SQIQL_REF_IMPL)_testvectors: $(SQIQL_REF_HOST_AUDIT_LIB)
