# Makefile for tgraph-gtk
# simple installation for a small Python project

# by default use ~/.local for normal users, /usr/local as root
PREFIX ?= $(shell if [ "$$(id -u)" = 0 ]; then echo "/usr/local"; else echo "$$HOME/.local"; fi)
BINDIR  = $(PREFIX)/bin

SRC_DIR = src

PYTHON = python3

# list of files we install
SCRIPT_SRC = $(SRC_DIR)/tgraph.py
MODULE_SRC = $(SRC_DIR)/tdata.py
TEXT_SRC   = $(SRC_DIR)/tgraph.txt

.PHONY: all install uninstall clean help

all:  ## default target does nothing
	@echo "run 'make install' to install tgraph."

install:  ## install into $(PREFIX)
	@echo "Installing to prefix $(PREFIX)"
	install -d "$(BINDIR)"
	install -m 755 "$(SCRIPT_SRC)" "$(BINDIR)/tgraph"
	install -m 644 "$(MODULE_SRC)" "$(BINDIR)/tdata.py"
	# documentation goes under share/doc/tgraph
	install -d "$(PREFIX)/share/doc/tgraph"
	install -m 644 "$(TEXT_SRC)" "$(PREFIX)/share/doc/tgraph/tgraph.txt"

uninstall:  ## remove installed files
	@echo "Removing files from $(BINDIR) and documentation"
	rm -f "$(BINDIR)/tgraph" "$(BINDIR)/tdata.py"
	rm -f "$(PREFIX)/share/doc/tgraph/tgraph.txt"
	rmdir --ignore-fail-on-non-empty "$(PREFIX)/share/doc/tgraph"

clean:  ## nothing to clean for now
	@echo "nothing to clean"

check:  ## run the test suite (requires pytest)
	@$(PYTHON) -m pytest tests

help:  ## show this help
	@grep -E '^[a-zA-Z_-]+:.*?##' $(MAKEFILE_LIST) | \
	sort | awk 'BEGIN {FS=":"} {printf "%-20s %s\n", $$1, $$2}'
