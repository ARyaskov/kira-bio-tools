import argparse
import difflib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def normalized(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    prefixes = ("##bcftools", "##kira_bt", "##fileDate=", "# This file was produced by", "# The command line was:")
    return [line for line in lines if not line.startswith(prefixes)]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--kira", type=Path, required=True)
    parser.add_argument("--bcftools", type=Path, required=True)
    parser.add_argument("--samtools", type=Path, required=True)
    parser.add_argument("--shell", type=Path, default=Path("bash"))
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--modes", nargs="*")
    parser.add_argument("--cases", nargs="*")
    parser.add_argument("--timeout", type=int, default=60)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="bcftools-compat-", dir=args.work_dir)).resolve()
    source = Path(__file__).resolve().parent
    fixtures = sorted(path for path in source.glob("*/test*") if path.is_dir() and
                      (not args.modes or path.parent.name in args.modes) and
                      (not args.cases or path.relative_to(source).as_posix() in {case.replace("\\", "/") for case in args.cases}))
    if not fixtures:
        raise RuntimeError("no fixtures selected")
    env = os.environ.copy()
    paths = [args.kira.resolve().parent, args.bcftools.resolve().parent,
             args.samtools.resolve().parent, args.shell.resolve().parent]
    env["PATH"] = os.pathsep.join(map(str, paths)) + os.pathsep + env.get("PATH", "")
    env["KIRA_BT_NO_PROGRESS"] = "1"
    env["MSYS2_ARG_CONV_EXCL"] = "*"
    versions = {}
    for name, executable in [("kira", args.kira), ("bcftools", args.bcftools), ("samtools", args.samtools)]:
        result = subprocess.run([str(executable), "--version"], env=env, capture_output=True, timeout=args.timeout)
        versions[name] = {"path": str(executable.resolve()), "exit_code": result.returncode,
                          "stdout": result.stdout.decode(encoding="utf-8", errors="replace")}
        if result.returncode:
            raise RuntimeError(f"{name} version check failed")
    (work / "versions.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")
    receipts = []
    for number, fixture in enumerate(fixtures, 1):
        relative = fixture.relative_to(source)
        target = work / relative
        shutil.copytree(fixture, target, ignore=shutil.ignore_patterns("*.ani", "*.ani.meta", "*.ktile"))
        if not target.resolve().is_relative_to(work):
            raise RuntimeError("fixture leaves isolated workspace")
        for path in target.iterdir():
            if path.is_symlink():
                raise RuntimeError(f"symlink in fixture: {relative}")
            if path.is_file() and (".kira.vcf" in path.name or ".bcf.vcf" in path.name):
                path.unlink()
            elif path.suffix == ".sh":
                path.write_text(path.read_text(encoding="utf-8"), encoding="utf-8", newline="\n")
        row = {"fixture": str(relative), "status": "NOT_RUN", "commands": []}
        for script in ["bcftools.sh", "kira.sh"]:
            if not (target / script).exists():
                row["status"] = "MISSING_SCRIPT"
                break
            try:
                result = subprocess.run([str(args.shell), "--noprofile", "--norc", "-eo", "pipefail", script],
                                        cwd=target, env=env, capture_output=True, timeout=args.timeout)
                (target / (script + ".stdout.log")).write_bytes(result.stdout)
                (target / (script + ".stderr.log")).write_bytes(result.stderr)
                row["commands"].append({"script": script, "exit_code": result.returncode})
                if result.returncode:
                    row["status"] = "EXECUTION_FAILED"
                    break
            except subprocess.TimeoutExpired:
                row["status"] = "TIMEOUT"
                break
        else:
            actual = sorted(target.glob("*.kira.vcf"))
            expected = sorted(target.glob("*.bcf.vcf"))
            if len(actual) != 1 or len(expected) != 1:
                row["status"] = "MISSING_OUTPUT"
            else:
                a, b = normalized(expected[0]), normalized(actual[0])
                row["status"] = "PASS" if a == b else "DIFFERENT"
                row["data_equal"] = [line for line in a if not line.startswith("#")] == [line for line in b if not line.startswith("#")]
                if a != b:
                    (target / "oracle.diff").write_text("\n".join(difflib.unified_diff(a, b, fromfile=expected[0].name, tofile=actual[0].name)), encoding="utf-8")
        receipts.append(row)
        (work / "results.json").write_text(json.dumps(receipts, indent=2), encoding="utf-8")
        if number % 10 == 0 or row["status"] != "PASS" or number == len(fixtures):
            print(f"{number}/{len(fixtures)} {relative}: {row['status']}", flush=True)
    summary = {status: sum(row["status"] == status for row in receipts) for status in sorted({row["status"] for row in receipts})}
    print(json.dumps({"work_dir": str(work), "summary": summary}), flush=True)
    return 0 if all(row["status"] == "PASS" for row in receipts) else 1


if __name__ == "__main__":
    raise SystemExit(main())
