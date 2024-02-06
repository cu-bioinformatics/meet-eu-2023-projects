#!/usr/bin/env python
import argparse
import pathlib
import shutil
import subprocess
import sys
from typing import List

from utils.utils import red_text

SCRIPT_NAME = pathlib.Path(__file__).name

DEPENDANCIES = ["parallel"]
PREFIX = "predicted_"
PREDICT_SCRIPT_NAME = "predict.py"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "target_features", type=str, help="path to csv file with target features"
    )
    parser.add_argument("model", type=str, help="path to pickled model")
    parser.add_argument(
        "files",
        type=str,
        nargs="*",
        help="files with ligand features on which to run parallel predictions.",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        default=PREFIX,
        help="prefix that will be added to each resulting temporary file (default: 'predicted_')",
    )
    parser.add_argument(
        "-p", "--processes", type=int, default=1, help="number of threads to use."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default=".",
        help="directory to which outputs should be written.",
    )
    return parser.parse_args()


def missing_dependancies() -> List[str]:
    missing = []
    for dependancy in DEPENDANCIES:
        if shutil.which(dependancy) is None:
            missing.append(dependancy)
    return missing


def main() -> int:
    args = parse_args()

    # Checking the paths chore.
    missing = missing_dependancies()
    if missing:
        print(
            f"{SCRIPT_NAME}: {red_text('error')}: missing dependancies: ({', '.join(missing)})"
        )
        return 2

    predict_script = pathlib.Path(__file__).resolve().parent / PREDICT_SCRIPT_NAME
    predict_script = predict_script
    if not predict_script.exists():
        print(
            f"{SCRIPT_NAME}: {red_text('error')}: {PREDICT_SCRIPT_NAME} not found ({predict_script})"
        )
        return 1

    target_features_path = pathlib.Path(args.target_features)
    model_path = pathlib.Path(args.model)
    if not target_features_path.exists():
        print(
            f"{red_text(SCRIPT_NAME)}: {red_text('error')}: {target_features_path} No such file or directory!"
        )
        return 1
    if not model_path.exists():
        print(
            f"{red_text(SCRIPT_NAME)}: {red_text('error')}: {model_path} No such file or directory!"
        )
        return 1

    files = [pathlib.Path(file) for file in args.files]
    files_doesnt_exist = [not file.exists() for file in files]
    if any(files_doesnt_exist):
        print(
            f"{SCRIPT_NAME}: {red_text('error')}: {files[files_doesnt_exist.index(True)]}: No such file or directory!"
        )
        return 1

    outdir = pathlib.Path(args.outdir)
    if outdir.exists() and not outdir.is_dir():
        print(f"{SCRIPT_NAME}: {red_text('error')}: {outdir} is not a directory!")
        return 1
    if not outdir.exists():
        outdir.mkdir(parents=True)

    # The actual scripts body.
    escape = lambda p: p.as_posix().replace(" ", r"\ ")
    command = (
        f'parallel -j {args.processes} "{escape(predict_script)}" '
        f'{{}} "{escape(target_features_path)}" "{escape(model_path)}" '
        f'"{outdir}/{args.prefix}"{{/}} ::: {" ".join(map(escape, files))}'
    )

    return subprocess.run(command, shell=True).returncode


if __name__ == "__main__":
    sys.exit(main())
