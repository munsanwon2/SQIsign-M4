#!/usr/bin/env python3
"""Fail-closed audit for dependencies unsuitable for a bare-metal M4 build.

Legacy imports are scanned conservatively as text. Curated imports can pass a
library manifest; in that mode every manifested production translation unit is
preprocessed with its TARGET_ARM definitions and only active backend code is
scanned. This excludes deliberately hosted ``#else`` branches without hiding
forbidden tokens that are active on Cortex-M4.
"""

from __future__ import annotations

from pathlib import Path, PurePosixPath
import argparse
import json
import re
import subprocess
import sys


SUFFIXES = {".c", ".h", ".S", ".s", ".inc", ".in", ".inl"}
CODE_RULES = {
    "GMP/mini-GMP": r"\b(?:mpz|mpq|mpf|__gmp)[A-Za-z0-9_]*",
    "host stdio": r"\b(?:fopen|freopen|fread|fwrite|fprintf|printf|puts|perror|FILE)\b",
    "host process/environment": r"\b(?:system|popen|fork|exec[lvpe]*|getenv|setenv|unsetenv|exit|abort)\b",
    "host heap allocation": r"\b(?:malloc|calloc|realloc|free)\s*\(",
    "host heap syscall": r"\b_sbrk\b",
    "unsupported allocator API": r"\b(?:aligned_alloc|posix_memalign|memalign|valloc|pvalloc)\s*\(",
    "stack allocation": r"\b(?:alloca|__builtin_alloca)\s*\(",
    "threads": r"\b(?:pthread_|thrd_|mtx_)",
    "dynamic loader": r"\b(?:dlopen|dlsym|dlclose)\b",
}
INCLUDE_RULES = {
    "GMP/mini-GMP": {"gmp.h", "mini-gmp.h"},
    "host stdio": {"stdio.h"},
    "threads": {"pthread.h", "threads.h"},
    "dynamic loader": {"dlfcn.h"},
}
INCLUDE_RE = re.compile(r"^\s*#\s*include\s*[<\"]([^>\"]+)[>\"]", re.M)
LINE_MARKER_RE = re.compile(r'^\s*#\s+\d+\s+"([^"]+)"')
RANDOMBYTES_DEF = re.compile(
    r"\b(?:int|void|size_t)\s+randombytes(?:_[A-Za-z0-9_]+)?\s*\([^;{}]*\)\s*\{",
    re.S,
)


class AuditError(RuntimeError):
    pass


def strip_comments(source: str) -> str:
    """Remove C comments while preserving newlines and string literals."""
    out: list[str] = []
    index = 0
    state = "code"
    while index < len(source):
        char = source[index]
        next_char = source[index + 1] if index + 1 < len(source) else ""
        if state == "code":
            if char == "/" and next_char == "*":
                state = "block"
                out.extend("  ")
                index += 2
                continue
            if char == "/" and next_char == "/":
                state = "line"
                out.extend("  ")
                index += 2
                continue
            if char == '"':
                state = "string"
            elif char == "'":
                state = "char"
            out.append(char)
        elif state == "block":
            out.append("\n" if char == "\n" else " ")
            if char == "*" and next_char == "/":
                out.append(" ")
                index += 2
                state = "code"
                continue
        elif state == "line":
            out.append("\n" if char == "\n" else " ")
            if char == "\n":
                state = "code"
        else:
            out.append(char)
            if char == "\\" and next_char:
                out.append(next_char)
                index += 2
                continue
            if (state == "string" and char == '"') or (state == "char" and char == "'"):
                state = "code"
        index += 1
    return "".join(out)


def strip_literals(source: str) -> str:
    """Remove C string/character literal contents after comments are gone."""
    out: list[str] = []
    index = 0
    quote: str | None = None
    while index < len(source):
        char = source[index]
        if quote is None:
            if char in {'"', "'"}:
                quote = char
                out.append(" ")
            else:
                out.append(char)
        else:
            out.append("\n" if char == "\n" else " ")
            if char == "\\" and index + 1 < len(source):
                index += 1
                out.append("\n" if source[index] == "\n" else " ")
            elif char == quote:
                quote = None
        index += 1
    return "".join(out)


def scan_active_text(source: str, label: str, included_headers: set[str] | None = None) -> list[str]:
    errors: list[str] = []
    without_comments = strip_comments(source)
    for header in INCLUDE_RE.findall(without_comments):
        basename = PurePosixPath(header).name
        for rule, forbidden in INCLUDE_RULES.items():
            if basename in forbidden:
                errors.append(f"{label}: {rule}")
    for basename in included_headers or set():
        for rule, forbidden in INCLUDE_RULES.items():
            if basename in forbidden:
                errors.append(f"{label}: {rule}")
    code = strip_literals(without_comments)
    if "/dev/urandom" in without_comments or re.search(r"\bgetrandom\b", code):
        errors.append(f"{label}: host entropy")
    for rule, pattern in CODE_RULES.items():
        if re.search(pattern, code):
            errors.append(f"{label}: {rule}")
    if RANDOMBYTES_DEF.search(code):
        errors.append(f"{label}: bundled randombytes function definition")
    return errors


def is_within(path: Path, root: Path) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def backend_preprocessed_text(output: str, root: Path, cwd: Path) -> tuple[str, set[str]]:
    active: list[str] = []
    included_headers: set[str] = set()
    in_backend = False
    for line in output.splitlines(keepends=True):
        marker = LINE_MARKER_RE.match(line)
        if marker:
            name = marker.group(1)
            if name.startswith("<") and name.endswith(">"):
                in_backend = False
                continue
            candidate = Path(name)
            if not candidate.is_absolute():
                candidate = cwd / candidate
            try:
                candidate = candidate.resolve()
            except OSError:
                in_backend = False
                continue
            in_backend = is_within(candidate, root)
            if not in_backend:
                included_headers.add(candidate.name)
            continue
        if in_backend:
            active.append(line)
    return "".join(active), included_headers


def manifest_values(manifest: dict[str, object], name: str) -> list[str]:
    value = manifest.get(name)
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise AuditError(f"manifest {name} must be a string list")
    return value


def audit_curated(root: Path, manifest_path: Path, compiler: str) -> list[str]:
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise AuditError(f"cannot read curated manifest: {manifest_path}: {exc}") from exc
    if not isinstance(manifest, dict):
        raise AuditError("curated manifest root must be an object")
    definitions = manifest_values(manifest, "pqm4_compile_definitions")
    if "TARGET_ARM" not in {item.split("=", 1)[0] for item in definitions}:
        raise AuditError("curated firmware audit requires TARGET_ARM in manifest definitions")
    include_directories = manifest_values(manifest, "pqm4_compile_include_directories")
    units = manifest.get("translation_units")
    if not isinstance(units, list) or not units:
        raise AuditError("manifest translation_units must be a non-empty list")

    base_command = [
        compiler,
        "-E",
        "-mcpu=cortex-m4",
        "-mthumb",
        "-ffreestanding",
        "-std=c11",
    ]
    base_command.extend(f"-D{definition}" for definition in definitions)
    for relative in include_directories:
        path = PurePosixPath(relative)
        if path.is_absolute() or ".." in path.parts:
            raise AuditError(f"unsafe manifested include directory: {relative!r}")
        include = (root / Path(*path.parts)).resolve()
        if not is_within(include, root) or not include.is_dir():
            raise AuditError(f"missing manifested include directory: {relative}")
        base_command.append(f"-I{include}")

    errors: list[str] = []
    for record in units:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise AuditError("invalid translation_units record")
        relative = PurePosixPath(record["path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise AuditError(f"unsafe manifested translation unit: {relative}")
        source = (root / Path(*relative.parts)).resolve()
        if not is_within(source, root) or not source.is_file():
            raise AuditError(f"missing manifested translation unit: {relative}")
        language = "assembler-with-cpp" if source.suffix == ".S" else "c"
        command = [*base_command, "-x", language, str(source)]
        completed = subprocess.run(
            command,
            cwd=root,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding="utf-8",
            errors="replace",
            check=False,
        )
        if completed.returncode != 0:
            detail = completed.stderr.strip().splitlines()
            message = detail[-1] if detail else f"compiler exited {completed.returncode}"
            raise AuditError(f"preprocessing failed for {relative}: {message}")
        active, included = backend_preprocessed_text(completed.stdout, root, root)
        errors.extend(scan_active_text(active, relative.as_posix(), included))
    return errors


def audit_legacy(root: Path) -> list[str]:
    errors: list[str] = []
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix not in SUFFIXES:
            continue
        source = path.read_text(errors="replace")
        path_errors = scan_active_text(source, str(path))
        if path.name == "malloc_compat.c":
            # The overlay deliberately exports the four libc allocator names
            # so imported code is redirected into the bounded qlw arena.  It
            # is not a hosted heap dependency, but accept it only in this
            # mechanically narrow form; a stray/native allocator call must
            # still fail the conservative source audit.
            code = strip_literals(strip_comments(source))
            exact_routes = {
                "malloc": r"void\s*\*\s*malloc\s*\(\s*size_t\s+size\s*\)\s*\{\s*return\s+qlw_malloc\s*\(\s*size\s*\)\s*;\s*\}",
                "calloc": r"void\s*\*\s*calloc\s*\(\s*size_t\s+count\s*,\s*size_t\s+size\s*\)\s*\{\s*return\s+qlw_calloc\s*\(\s*count\s*,\s*size\s*\)\s*;\s*\}",
                "realloc": r"void\s*\*\s*realloc\s*\(\s*void\s*\*\s*ptr\s*,\s*size_t\s+size\s*\)\s*\{\s*return\s+qlw_realloc\s*\(\s*ptr\s*,\s*size\s*\)\s*;\s*\}",
                "free": r"void\s+free\s*\(\s*void\s*\*\s*ptr\s*\)\s*\{\s*qlw_free\s*\(\s*ptr\s*\)\s*;\s*\}",
            }
            valid = (
                '"qlapoti_workspace.h"' in source
                and "<stdlib.h>" not in source
                and all(len(re.findall(pattern, code, re.S)) == 1
                        for pattern in exact_routes.values())
                and len(re.findall(r"\b(?:malloc|calloc|realloc|free)\s*\(", code)) == 4
            )
            if valid:
                path_errors = [
                    error for error in path_errors
                    if not error.endswith(": host heap allocation")
                ]
            else:
                path_errors.append(f"{path}: invalid bounded allocator shim")
        errors.extend(path_errors)
        if path.name.lower() in {
            "randombytes.c",
            "randombytes_system.c",
            "randombytes.s",
        }:
            errors.append(f"{path}: bundled randombytes implementation file")
    return errors


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--library-manifest", type=Path)
    parser.add_argument("--cc", default="arm-none-eabi-gcc")
    args = parser.parse_args()
    root = args.root.resolve()
    if not root.is_dir():
        print(f"ERROR: not a directory: {root}", file=sys.stderr)
        return 2
    try:
        if args.library_manifest is not None:
            errors = audit_curated(root, args.library_manifest.resolve(), args.cc)
        else:
            errors = audit_legacy(root)
    except (AuditError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if errors:
        for error in sorted(set(errors)):
            print("ERROR:", error, file=sys.stderr)
        return 1
    mode = "ACTIVE TARGET_ARM" if args.library_manifest is not None else "CONSERVATIVE SOURCE"
    print(f"NO OBVIOUS HOST-ONLY/FOREIGN-RNG DEPENDENCIES ({mode}): PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
