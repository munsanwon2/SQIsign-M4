#!/usr/bin/env python3
"""Create deterministic pqm4-oriented manifests from SQIsign's CMake graph.

The tool asks CMake's File API for the transitive source closure of the real
``sqisign_lvl1``, ``sqisign_lvl3``, or ``sqisign_lvl5`` production library.
It does not infer the library from directory names.  A small, explicit pqm4
policy then removes implementations that must be supplied by the platform
(host RNG and host heap/timing helpers), rejects tests/NIST/GMP contamination,
and optionally exports a directory-preserving source tree.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile


SCHEMA_VERSION = 2
TRANSLATION_UNIT_SUFFIXES = {".c", ".s", ".S"}
HEADER_SUFFIXES = {".h"}
LEVELS = (1, 3, 5)

PQM4_EXCLUSIONS = {
    "src/common/generic/randombytes_system.c": (
        "host OS RNG; pqm4 must provide its hardware/platform randombytes adapter"
    ),
    "src/common/generic/tools.c": (
        "host clock/stdio benchmarking helpers; no production library caller requires them"
    ),
    "src/common/generic/mem.c": (
        "host free()-based secure cleanup; pqm4 must provide the no-heap secure-clear adapter"
    ),
    "src/quaternion/ref/generic/test/random_input_generation.c": (
        "test-only random fixture with malloc/stdio; currently listed by the upstream production CMake target but has no production caller"
    ),
}

REQUIRED_QLAPOTI_SOURCES = {
    "src/id2iso/ref/lvlx/dim2id2iso.c",
    "src/quaternion/ref/generic/dim2.c",
    "src/quaternion/ref/generic/integers.c",
    "src/quaternion/ref/generic/qlapoti.c",
    "src/quaternion/ref/generic/rationals.c",
    "src/signature/ref/lvlx/keygen.c",
    "src/signature/ref/lvlx/sign.c",
}

FORBIDDEN_PATH_COMPONENTS = {
    "apps",
    "bench-results",
    "kat",
    "mini-gmp",
    "nistapi",
    "test",
    "tests",
}

GMP_INCLUDE = re.compile(
    r"^\s*#\s*include\s*[<\"][^\">]*gmp[^\">]*[\">]", re.IGNORECASE | re.MULTILINE
)
INCLUDE_DIRECTIVE = re.compile(
    r"^\s*#\s*include\s*([<\"])([^\">]+)[\">]", re.MULTILINE
)
DEFINE_FRAGMENT = re.compile(r"(?:^|\s)-D([A-Za-z_][A-Za-z0-9_]*(?:=[^\s]+)?)")

# CMake's Release configuration also contributes optimization, diagnostics,
# and host-build policy.  pqm4 owns optimization and warning policy, so only
# target-independent options with source-level semantic relevance cross the
# export boundary.  The complete CMake option set remains recorded for audit.
PQM4_COMPILE_OPTION_ALLOWLIST = {
    "-fno-strict-aliasing",
    "-fvisibility=hidden",
    "-funroll-loops",
}


class ManifestError(RuntimeError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_relative(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError as exc:
        raise ManifestError(f"CMake selected path outside source tree: {path}") from exc


def resolve_cmake_path(raw_path: str, source_root: Path) -> Path:
    path = Path(raw_path)
    if not path.is_absolute():
        path = source_root / path
    return path.resolve()


def forbidden_path(path: str) -> bool:
    return any(part.lower() in FORBIDDEN_PATH_COMPONENTS for part in Path(path).parts)


def file_record(path: Path, source_root: Path, **extra: object) -> dict[str, object]:
    content = path.read_bytes()
    record: dict[str, object] = {
        "path": normalized_relative(path, source_root),
        "bytes": len(content),
        "sha256": hashlib.sha256(content).hexdigest(),
    }
    record.update(extra)
    return record


def configure_codemodel(source_root: Path, cmake: str) -> tuple[Path, dict[str, object], tempfile.TemporaryDirectory[str]]:
    temporary = tempfile.TemporaryDirectory(prefix="sqisign-pqm4-codemodel-")
    build_root = Path(temporary.name)
    query = build_root / ".cmake/api/v1/query/codemodel-v2"
    query.parent.mkdir(parents=True)
    query.touch()

    command = [
        cmake,
        "-S", str(source_root),
        "-B", str(build_root),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_SYSTEM_NAME=Generic",
        "-DCMAKE_SYSTEM_PROCESSOR=arm",
        "-DCMAKE_TRY_COMPILE_TARGET_TYPE=STATIC_LIBRARY",
        "-DENABLE_TESTS=OFF",
        "-DENABLE_KAT_TESTS=OFF",
        "-DENABLE_SIGN=ON",
        "-DENABLE_GMP=OFF",
        "-DSQISIGN_BUILD_TYPE=ref",
        "-DGF_RADIX=32",
    ]
    process = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if process.returncode != 0:
        temporary.cleanup()
        raise ManifestError(
            "CMake codemodel configuration failed\n"
            f"command: {' '.join(command)}\n"
            f"stdout:\n{process.stdout}\nstderr:\n{process.stderr}"
        )

    reply_root = build_root / ".cmake/api/v1/reply"
    indexes = sorted(reply_root.glob("index-*.json"))
    if len(indexes) != 1:
        temporary.cleanup()
        raise ManifestError(f"expected one CMake File API index, found {len(indexes)}")
    index = json.loads(indexes[0].read_text(encoding="utf-8"))
    try:
        codemodel_name = index["reply"]["codemodel-v2"]["jsonFile"]
    except (KeyError, TypeError) as exc:
        temporary.cleanup()
        raise ManifestError("CMake File API did not return codemodel-v2") from exc
    codemodel = json.loads((reply_root / codemodel_name).read_text(encoding="utf-8"))
    return reply_root, codemodel, temporary


def select_configuration(codemodel: dict[str, object]) -> dict[str, object]:
    configurations = codemodel.get("configurations", [])
    if not isinstance(configurations, list) or not configurations:
        raise ManifestError("CMake codemodel has no configurations")
    for configuration in configurations:
        if configuration.get("name") == "Release":
            return configuration
    if len(configurations) == 1:
        return configurations[0]
    raise ManifestError("CMake codemodel has no unambiguous Release configuration")


def load_target_graph(reply_root: Path, configuration: dict[str, object]) -> tuple[dict[str, dict[str, object]], dict[str, str]]:
    targets_by_id: dict[str, dict[str, object]] = {}
    ids_by_name: dict[str, str] = {}
    for reference in configuration.get("targets", []):
        target = json.loads((reply_root / reference["jsonFile"]).read_text(encoding="utf-8"))
        target_id = target["id"]
        target_name = target["name"]
        targets_by_id[target_id] = target
        if target_name in ids_by_name:
            raise ManifestError(f"duplicate CMake target name {target_name}")
        ids_by_name[target_name] = target_id
    return targets_by_id, ids_by_name


def target_closure(root_id: str, targets_by_id: dict[str, dict[str, object]]) -> list[dict[str, object]]:
    pending = [root_id]
    seen: set[str] = set()
    closure: list[dict[str, object]] = []
    while pending:
        target_id = pending.pop()
        if target_id in seen:
            continue
        if target_id not in targets_by_id:
            raise ManifestError(f"CMake dependency target is missing from codemodel: {target_id}")
        seen.add(target_id)
        target = targets_by_id[target_id]
        closure.append(target)
        pending.extend(dependency["id"] for dependency in target.get("dependencies", []))
    return sorted(closure, key=lambda target: target["name"])


def collect_headers(
    translation_units: list[Path],
    include_directories: list[Path],
    source_root: Path,
) -> list[Path]:
    """Return only project headers transitively reached by production TUs.

    Copying every file below a CMake include directory accidentally exports
    benchmark and test helper headers that are not part of the library.  This
    follows C/C++ include directives instead.  System headers are deliberately
    left to the pqm4 toolchain; only includes resolving inside ``source_root``
    enter the content-addressed export.
    """

    headers: set[Path] = set()
    pending: list[Path] = list(translation_units)
    scanned: set[Path] = set()
    while pending:
        current = pending.pop()
        current = current.resolve()
        if current in scanned:
            continue
        scanned.add(current)
        if current.suffix not in TRANSLATION_UNIT_SUFFIXES | HEADER_SUFFIXES:
            continue
        text = current.read_text(encoding="utf-8", errors="strict")
        for match in INCLUDE_DIRECTIVE.finditer(text):
            opening, include_name = match.groups()
            candidates: list[Path] = []
            if opening == '"':
                candidates.append(current.parent / include_name)
            candidates.extend(directory / include_name for directory in include_directories)
            resolved = next(
                (candidate.resolve() for candidate in candidates if candidate.is_file()),
                None,
            )
            if resolved is None:
                continue
            relative = normalized_relative(resolved, source_root)
            if resolved.suffix not in HEADER_SUFFIXES:
                raise ManifestError(
                    f"project include does not resolve to a header: {current}: {include_name}"
                )
            if forbidden_path(relative):
                raise ManifestError(
                    f"production header closure reaches a forbidden path: {relative}"
                )
            if resolved not in headers:
                headers.add(resolved)
                pending.append(resolved)
    return sorted(headers, key=lambda path: normalized_relative(path, source_root))


def manifest_for_level(
    level: int,
    source_root: Path,
    targets_by_id: dict[str, dict[str, object]],
    ids_by_name: dict[str, str],
) -> dict[str, object]:
    target_name = f"sqisign_lvl{level}"
    if target_name not in ids_by_name:
        raise ManifestError(f"CMake production target is missing: {target_name}")
    closure = target_closure(ids_by_name[target_name], targets_by_id)

    sources_by_path: dict[Path, set[str]] = defaultdict(set)
    include_directories: set[Path] = set()
    include_directory_order: list[Path] = []
    definitions: set[str] = set()
    compile_options: set[str] = set()
    target_records: list[dict[str, str]] = []

    for target in closure:
        target_records.append({"name": target["name"], "type": target["type"]})
        for compile_group in target.get("compileGroups", []):
            definitions.update(item["define"] for item in compile_group.get("defines", []))
            for fragment in compile_group.get("compileCommandFragments", []):
                fragment_text = fragment.get("fragment", "")
                definitions.update(DEFINE_FRAGMENT.findall(fragment_text))
                try:
                    fragment_tokens = shlex.split(fragment_text, posix=True)
                except ValueError as exc:
                    raise ManifestError(
                        f"could not tokenize CMake compile fragment for {target['name']}: "
                        f"{fragment_text!r}"
                    ) from exc
                compile_options.update(
                    token
                    for token in fragment_tokens
                    if not token.startswith("-D") and not token.startswith("-I")
                )
            for include in compile_group.get("includes", []):
                include_path = resolve_cmake_path(include["path"], source_root)
                normalized_relative(include_path, source_root)
                include_directories.add(include_path)
                if include_path not in include_directory_order:
                    include_directory_order.append(include_path)

        for source in target.get("sources", []):
            path = resolve_cmake_path(source["path"], source_root)
            if path.suffix not in TRANSLATION_UNIT_SUFFIXES:
                continue
            if not path.is_file():
                raise ManifestError(f"CMake-selected translation unit is missing: {path}")
            normalized_relative(path, source_root)
            sources_by_path[path].add(target["name"])

    raw_paths = {normalized_relative(path, source_root) for path in sources_by_path}
    missing_qlapoti = sorted(REQUIRED_QLAPOTI_SOURCES - raw_paths)
    required_level_sources = {
        f"src/precomp/ref/lvl{level}/endomorphism_action.c",
        f"src/precomp/ref/lvl{level}/quaternion_data.c",
    }
    missing_level = sorted(required_level_sources - raw_paths)
    if missing_qlapoti or missing_level:
        raise ManifestError(
            f"{target_name} lacks required production sources: "
            f"{missing_qlapoti + missing_level}"
        )

    contaminated = sorted(
        path
        for path in raw_paths
        if forbidden_path(path) and path not in PQM4_EXCLUSIONS
    )
    if contaminated:
        raise ManifestError(
            f"test/NIST/GMP/app source entered {target_name} closure: {contaminated}"
        )

    selected_records: list[dict[str, object]] = []
    excluded_records: list[dict[str, object]] = []
    for path in sorted(sources_by_path, key=lambda item: normalized_relative(item, source_root)):
        relative = normalized_relative(path, source_root)
        targets = sorted(sources_by_path[path])
        if relative in PQM4_EXCLUSIONS:
            excluded_records.append(
                file_record(
                    path,
                    source_root,
                    reason=PQM4_EXCLUSIONS[relative],
                    selected_by=targets,
                )
            )
            continue
        text = path.read_text(encoding="utf-8", errors="strict") if path.suffix == ".c" else ""
        if GMP_INCLUDE.search(text) or "GMP_LIMB_BITS" in text:
            raise ManifestError(f"GMP-dependent source entered fixed manifest: {relative}")
        selected_records.append(
            file_record(path, source_root, selected_by=targets)
        )

    excluded_paths = {record["path"] for record in excluded_records}
    if "src/common/generic/randombytes_system.c" not in excluded_paths:
        raise ManifestError("CMake graph no longer exposes the host RNG exclusion explicitly")

    existing_include_directories = {
        path for path in include_directories if path.is_dir()
    }
    missing_include_directories = sorted(
        normalized_relative(path, source_root)
        for path in include_directories
        if not path.is_dir()
    )
    selected_paths = [source_root / record["path"] for record in selected_records]
    headers = collect_headers(
        selected_paths,
        [path for path in include_directory_order if path in existing_include_directories],
        source_root,
    )
    for header in headers:
        header_text = header.read_text(encoding="utf-8", errors="strict")
        if GMP_INCLUDE.search(header_text):
            raise ManifestError(
                f"GMP-dependent header entered fixed manifest: "
                f"{normalized_relative(header, source_root)}"
            )
    header_records = [file_record(path, source_root) for path in headers]
    used_include_directories = {
        directory
        for directory in existing_include_directories
        if any(
            header == directory or directory in header.parents
            for header in headers
        )
    }
    include_records = sorted(
        normalized_relative(path, source_root) for path in used_include_directories
    )
    pqm4_definitions = sorted(
        definition
        for definition in definitions
        if definition != "RANDOMBYTES_SYSTEM"
        and not definition.startswith("SQISIGN_TEST_REPS=")
    )
    pqm4_compile_options = sorted(
        option for option in compile_options if option in PQM4_COMPILE_OPTION_ALLOWLIST
    )

    return {
        "schema_version": SCHEMA_VERSION,
        "kind": "sqisign-pqm4-library-source-manifest",
        "level": level,
        "root_target": target_name,
        "selection": {
            "mechanism": "CMake File API codemodel-v2 transitive target closure",
            "configuration": "Release",
            "cmake_options": [
                "CMAKE_SYSTEM_NAME=Generic",
                "CMAKE_SYSTEM_PROCESSOR=arm",
                "ENABLE_TESTS=OFF",
                "ENABLE_KAT_TESTS=OFF",
                "ENABLE_SIGN=ON",
                "ENABLE_GMP=OFF",
                "SQISIGN_BUILD_TYPE=ref",
                "GF_RADIX=32",
            ],
            "target_closure": target_records,
        },
        "translation_units": selected_records,
        "headers": header_records,
        "include_directories": include_records,
        "cmake_missing_include_directories": missing_include_directories,
        "cmake_compile_definitions": sorted(definitions),
        "pqm4_compile_definitions": pqm4_definitions,
        "cmake_compile_options": sorted(compile_options),
        "pqm4_compile_options": pqm4_compile_options,
        "pqm4_compile_include_directories": include_records,
        "pqm4_compile_policy": {
            "definitions": (
                "CMake target-closure definitions except host RNG and test-count macros"
            ),
            "options": (
                "only target-independent semantic options; pqm4 retains ownership of "
                "optimization, language-standard, diagnostics, ABI, and CPU flags"
            ),
            "include_directories": (
                "existing source-tree include directories selected by the CMake target closure"
            ),
        },
        "excluded_translation_units": excluded_records,
        "external_platform_requirements": [
            {
                "interface": "sqisign_platform_randombytes(unsigned char *, size_t)",
                "provider": (
                    "pqm4 adapter with the exact size_t ABI; it may call pqm4 "
                    "randombytes(uint8_t *, size_t)"
                ),
            },
            {
                "interface": "sqisign_secure_clear/sqisign_secure_free",
                "provider": "pqm4 no-heap secure cleanup adapter if referenced",
            },
        ],
        "invariants": {
            "required_qlapoti_sources_present": True,
            "level_precomputation_present": True,
            "tests_and_benchmarks_absent": True,
            "nist_wrappers_absent": True,
            "gmp_sources_and_includes_absent": True,
            "host_rng_absent": True,
        },
    }


def write_output(output: Path, manifests: list[dict[str, object]], source_root: Path, export: bool) -> None:
    if output.exists():
        raise ManifestError(f"output already exists; refusing to overwrite: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=f".{output.name}.tmp-", dir=output.parent))
    try:
        for manifest in manifests:
            level_root = staging / f"lvl{manifest['level']}"
            level_root.mkdir(parents=True)
            manifest_path = level_root / "library-manifest.json"
            manifest_path.write_text(
                json.dumps(manifest, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
                newline="\n",
            )
            for category in ("translation_units", "headers"):
                for record in manifest[category]:
                    source = source_root / record["path"]
                    if source.stat().st_size != record["bytes"] or sha256_file(source) != record["sha256"]:
                        raise ManifestError(
                            f"source changed while manifest was generated: {record['path']}"
                        )
                    if export:
                        destination = level_root / "source" / record["path"]
                        destination.parent.mkdir(parents=True, exist_ok=True)
                        shutil.copyfile(source, destination)
                        if sha256_file(destination) != record["sha256"]:
                            raise ManifestError(
                                f"export hash mismatch after copy: {record['path']}"
                            )
        os.replace(staging, output)
    except Exception:
        shutil.rmtree(staging, ignore_errors=True)
        raise


def parse_levels(values: list[int]) -> list[int]:
    levels = sorted(set(values))
    unsupported = [level for level in levels if level not in LEVELS]
    if unsupported:
        raise ManifestError(f"unsupported levels: {unsupported}; expected 1, 3, or 5")
    return levels


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True, help="merged SQIsign source root")
    parser.add_argument("--output", type=Path, required=True, help="new manifest/export directory")
    parser.add_argument("--levels", type=int, nargs="+", default=list(LEVELS))
    parser.add_argument("--export", action="store_true", help="copy selected sources and headers")
    parser.add_argument("--cmake", default=os.environ.get("CMAKE", "cmake"))
    args = parser.parse_args()

    source_root = args.source.resolve()
    if not (source_root / "CMakeLists.txt").is_file():
        raise ManifestError(f"not an SQIsign CMake source tree: {source_root}")
    levels = parse_levels(args.levels)
    output = args.output.resolve()
    if output.exists():
        raise ManifestError(f"output already exists; refusing to overwrite: {output}")

    reply_root, codemodel, temporary = configure_codemodel(source_root, args.cmake)
    try:
        configuration = select_configuration(codemodel)
        targets_by_id, ids_by_name = load_target_graph(reply_root, configuration)
        manifests = [
            manifest_for_level(level, source_root, targets_by_id, ids_by_name)
            for level in levels
        ]
        write_output(output, manifests, source_root, args.export)
    finally:
        temporary.cleanup()

    for manifest in manifests:
        print(
            f"lvl{manifest['level']}: "
            f"translation_units={len(manifest['translation_units'])} "
            f"headers={len(manifest['headers'])} "
            f"excluded={len(manifest['excluded_translation_units'])}"
        )
    print("PQM4 LIBRARY SOURCE MANIFEST: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ManifestError as error:
        raise SystemExit(f"PQM4 LIBRARY SOURCE MANIFEST: FAIL\n{error}")
