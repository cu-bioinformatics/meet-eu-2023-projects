#!/usr/bin/env python
"""script for fetching data from ChEMBL."""
import argparse
import datetime
import pathlib
import sys
from typing import Any, List, Optional, Union

import requests
import pandas as pd

SCRIPT_NAME = pathlib.Path(__file__).name

BASE_URL = "https://www.ebi.ac.uk"
DATA_ENDPOINT = f"{BASE_URL}/chembl/api/data"

COMMANDS = ["assays", "compounds"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="A script for downloading bioactivity assay data from ChEMBL"
    )
    parser.add_argument(
        "command",
        type=str,
        help=f"command specyfing what to fetch (possible values: {', '.join(COMMANDS)})",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="path to a file, to which the output will be wriiten.",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="don't ask for confirmation when overwriting existing files.",
    )
    return parser.parse_args()


def fetch_url(endpoint: str, format: str = "json", **params) -> str:
    url = f"{DATA_ENDPOINT}/{endpoint}?format={format}"
    for name, value in params.items():
        url += f"&{name}={value}"
    return url


def next_page(response_body: dict[str, Any]) -> Optional[str]:
    return response_body["page_meta"]["next"]


def to_df(r: dict, offset: int) -> pd.DataFrame:
    content = {} 
    for mol in r["molecules"]:
        properties = mol["molecule_properties"]
        structures = mol["molecule_structures"]
        if properties is not None and structures is not None:
            content[mol["molecule_chembl_id"]] = {
                **properties,
                "canonical_smiles": structures["canonical_smiles"]
            }
    try:
        return pd.DataFrame(content).transpose()
    except ValueError:
        filename = f"{SCRIPT_NAME}_{datetime.datetime.now()}_unparsed.log"
        with open(filename, "a") as logfile:
            logfile.write(str({offset: content}))
    return pd.DataFrame()


def get_viruses_ids(batch_size: int = 500) -> List[Union[str, int]]:
    url = fetch_url("organism", limit=f"{batch_size}", l1="Viruses")
    ids = []
    response = requests.get(url)
    if response.ok:
        response = response.json()
        next = next_page(response)
        ids = [organism["tax_id"] for organism in response["organisms"]]
        while next is not None:
            next_page_body = requests.get(f"{BASE_URL}{next}").json()
            ids.extend([organism["tax_id"] for organism in next_page_body["organisms"]])
            next = next_page(next_page_body)
    return ids


def get_asssays() -> pd.DataFrame:
    url = fetch_url(
        "assay",
        limit=500,
        assay_type="B",
        description__icontains="Helicase",
        assay_organism__regex=".?virus.?",
    )
    response = requests.get(url)
    frames = []
    if response.ok:
        response_body = response.json()
        df = pd.DataFrame(
            response_body["assays"], columns=response_body["assays"][0].keys()
        )
        frames.append(df)
        next = next_page(response_body)
        # This does not work for now but that's not that importnat
        while next is not None:
            next_page_body = requests.get(f"{BASE_URL}{next}").json()
            frames.append(pd.DataFrame(next_page_body["assays"]))
            next = next_page(next_page_body)
    return pd.concat(frames, ignore_index=True)


def get_activities(assays_ids: List[str]):
    url = fetch_url("activity", limit=500, assay_chembl_id__in=",".join(assays_ids))
    response = requests.get(url)
    frames = []
    if response.ok:
        response_body = response.json()
        df = pd.DataFrame(
            response_body["activities"], columns=response_body["activities"][0].keys()
        )
        frames.append(df)
        next = next_page(response_body)
        # This does not work for now but that's not that importnat
        while next is not None:
            next_page_body = requests.get(f"{BASE_URL}{next}").json()
            frames.append(pd.DataFrame(next_page_body["activities"]))
            next = next_page(next_page_body)
    return pd.concat(frames, ignore_index=True)


def get_compounds(batch_size: int = 1000, output_filename: str = "ligands"):
    url = fetch_url("molecule", limit=1000)
    frames = []

    response = requests.get(url)
    response_body = response.json()
    frames.append(to_df(response_body, offset=0))
    with open(f"{output_filename}.tmp", "w") as file:
        file.write(frames[0].to_csv(header=True))

    next_page = response_body["page_meta"]["next"]
    total_batches = response_body["page_meta"]["total_count"] // batch_size + 1

    i = 1
    percent_done = i * 100 // total_batches
    while next_page is not None or i > total_batches:
        print(
            f"\033[2K\033[32m {'━' * int(percent_done * 0.8)}\033[31m{'━' * int((100 - percent_done) * 0.8)} "  # ]]]
            f"\033[0m{i}/{total_batches}",  # ]
            end="\r",
            flush=True,
        )
        response_body = requests.get(f"https://www.ebi.ac.uk{next_page}").json()
        next_page = response_body["page_meta"]["next"]
        offset = response_body["page_meta"]["offset"]
        frames.append(to_df(response_body, offset))
        with open(f"{output_filename}.tmp", "a") as file:
            file.write(frames[-1].to_csv(header=False))
        i += 1
        percent_done = i * 100 // total_batches
    return pd.concat(frames, ignore_index=True)


def fetch_assays(args: argparse.Namespace) -> int:
    assays = get_asssays()
    activities = get_activities(assays["assay_chembl_id"])

    if args.output:
        output_file = pathlib.Path(args.output)

        if output_file.is_dir():
            print(f"{red_text(SCRIPT_NAME)}: error: {output_file} is a directory.")
            return 1
        if (not args.force) and output_file.exists():
            print(f"[{yellow_text('warning')}] File {output_file} exists.", end=" ")
            choice = input("Do you want to overwrite it? [Y/n]: ")

            if choice.lower() not in ["y", "yes"]:
                print(activities.to_csv())
                return 0

        with open(output_file, "w") as file:
            file.write(activities.to_csv())
        return 0

    print(activities.to_csv())
    return 0


def fetch_compounds(args: argparse.Namespace) -> int:
    if args.output:
        compounds = get_compounds(output_filename=args.output)
    else:
        compounds = get_compounds()

    if args.output:
        output_file = pathlib.Path(args.output)

        if output_file.is_dir():
            print(f"{red_text(SCRIPT_NAME)}: error: {output_file} is a directory.")
            return 1
        if (not args.force) and output_file.exists():
            print(f"[{yellow_text('warning')}] File {output_file} exists.", end=" ")
            choice = input("Do you want to overwrite it? [Y/n]: ")

            if choice.lower() not in ["y", "yes"]:
                print(compounds.to_csv())
                return 0

        with open(output_file, "w") as file:
            file.write(compounds.to_csv())
        return 0

    print(compounds.to_csv())
    return 0


def main() -> int:
    args = parse_args()

    if args.command.lower() == "assays":
        return fetch_assays(args)

    if args.command.lower() == "compounds":
        return fetch_compounds(args)

    print(
        f"{red_text(SCRIPT_NAME)}: error: command must be one of: \033[1m{', '.join(COMMANDS)}\033[0m"
    )
    return 2


if __name__ == "__main__":
    sys.exit(main())
