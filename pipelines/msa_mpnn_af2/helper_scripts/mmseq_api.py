# Acknowledgment:
# The API code used in this project is adapted from the ColabFold repository by sokrypton.
# GitHub repository: https://github.com/sokrypton/ColabFold

import hashlib
import logging
import os
import random
import re
import tarfile
import time
from collections import Counter, OrderedDict
from typing import List, Tuple

import numpy as np
import requests
from tqdm import tqdm

logger = logging.getLogger(__name__)

TQDM_BAR_FORMAT = "{l_bar}{bar}| {n_fmt}/{total_fmt} [elapsed: {elapsed} remaining: {remaining}]"


def generate_fasta_sequence_from_pdb(pdb_file, output_folder=None):
    ca_pattern = re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")
    aa3to1 = {
        "ALA": "A",
        "VAL": "V",
        "PHE": "F",
        "PRO": "P",
        "MET": "M",
        "ILE": "I",
        "LEU": "L",
        "ASP": "D",
        "GLU": "E",
        "LYS": "K",
        "ARG": "R",
        "SER": "S",
        "THR": "T",
        "TYR": "Y",
        "HIS": "H",
        "CYS": "C",
        "ASN": "N",
        "GLN": "Q",
        "TRP": "W",
        "GLY": "G",
        "MSE": "M",
    }
    filename = os.path.basename(pdb_file).split(".")[0]
    chain_dict = dict()
    chain_list = []

    with open(pdb_file, "r") as fp:
        for line in fp:
            if line.startswith("ENDMDL"):
                break
            match_list = ca_pattern.findall(line)
            if match_list:
                resn = match_list[0][0] + match_list[0][2]
                chain = match_list[0][1] + match_list[0][3]
                if chain in chain_dict:
                    chain_dict[chain] += aa3to1[resn]
                else:
                    chain_dict[chain] = aa3to1[resn]
                    chain_list.append(chain)

    fasta_sequence = f">{filename}\n"
    for i, chain in enumerate(chain_list):
        fasta_sequence += chain_dict[chain]
        if i < len(chain_list) - 1:
            fasta_sequence += ":"

    if output_folder:
        output_file = os.path.join(output_folder, f"{filename}.fasta")
        with open(output_file, "w") as fp:
            fp.write(fasta_sequence)

    return chain_dict


def parse_fasta(fasta_string: str) -> Tuple[List[str], List[str]]:
    """Parses FASTA string and returns list of strings with amino-acid sequences.

    Arguments:
      fasta_string: The string contents of a FASTA file.

    Returns:
      A tuple of two lists:
      * A list of sequences.
      * A list of sequence descriptions taken from the comment lines. In the
        same order as the sequences.
    """
    sequences = []
    descriptions = []
    index = -1
    for line in fasta_string.splitlines():
        line = line.strip()
        if line.startswith("#"):
            continue
        if line.startswith(">"):
            index += 1
            descriptions.append(line[1:])  # Remove the '>' at the beginning.
            sequences.append("")
            continue
        elif not line:
            continue  # Skip blank lines.
        sequences[index] += line

    return sequences, descriptions


def get_hash(x):
    return hashlib.sha1(x.encode()).hexdigest()


def get_unique_sequences(seq_list):
    unique_seqs = list(OrderedDict.fromkeys(seq_list))
    return unique_seqs


def run_mmseq2_msa_rosettafold(seq, jobname, cov=50, id=90, max_msa=2048, mode="unpaired_paired", msa_name="msa"):
    assert mode in ["unpaired", "paired", "unpaired_paired"]
    seqs = [seq] if isinstance(seq, str) else seq

    # Process sequence
    seq = re.sub("[^A-Z:]", "", seq.replace("/", ":").upper())
    seq = re.sub(":+", ":", seq)
    seq = re.sub("^[:]+", "", seq)
    seq = re.sub("[:]+$", "", seq)
    seq = seq.replace(":", "/").split("/")

    print("Starting MSA...")
    get_msa(seq, jobname, cov=cov, id=id, max_msa=max_msa, mode=mode, msa_name=msa_name)

    print("Done!")


def get_msa(seq, jobname, cov=50, id=90, max_msa=2048, mode="unpaired_paired", msa_name="msa"):

    assert mode in ["unpaired", "paired", "unpaired_paired"]
    seqs = [seq] if isinstance(seq, str) else seq

    # collapse homooligomeric sequences
    counts = Counter(seqs)
    u_seqs = list(counts.keys())
    u_nums = list(counts.values())

    # expand homooligomeric sequences
    first_seq = "/".join(sum([[x] * n for x, n in zip(u_seqs, u_nums)], []))
    msa = [first_seq]

    path = os.path.join(jobname, f"{msa_name}")
    os.makedirs(path, exist_ok=True)
    if mode in ["paired", "unpaired_paired"] and len(u_seqs) > 1:
        print("getting paired MSA")
        out_paired = run_mmseqs2(u_seqs, f"{path}/", use_pairing=True)
        headers, sequences = [], []
        for a3m_lines in out_paired:
            n = -1
            for line in a3m_lines.split("\n"):
                if len(line) > 0:
                    if line.startswith(">"):
                        n += 1
                        if len(headers) < (n + 1):
                            headers.append([])
                            sequences.append([])
                        headers[n].append(line)
                    else:
                        sequences[n].append(line)
        # filter MSA
        with open(f"{path}/paired_in_{msa_name}.a3m", "w") as handle:
            for n, sequence in enumerate(sequences):
                handle.write(f">n{n}\n{''.join(sequence)}\n")
        os.system(f"hhfilter -i {path}/paired_in_{msa_name}.a3m -id {id} -cov {cov} -o {path}/paired_out_{msa_name}.a3m")
        with open(f"{path}/paired_out_{msa_name}.a3m", "r") as handle:
            for line in handle:
                if line.startswith(">"):
                    n = int(line[2:])
                    xs = sequences[n]
                    # expand homooligomeric sequences
                    xs = ["/".join([x] * num) for x, num in zip(xs, u_nums)]
                    msa.append("/".join(xs))

    if len(msa) < max_msa and (mode in ["unpaired", "unpaired_paired"] or len(u_seqs) == 1):
        print("getting unpaired MSA")
        out = run_mmseqs2(u_seqs, f"{path}/")
        Ls = [len(seq) for seq in u_seqs]
        sub_idx = []
        sub_msa = []
        sub_msa_num = 0
        for n, a3m_lines in enumerate(out):
            sub_msa.append([])
            with open(f"{path}/in_{msa_name}_{n}.a3m", "w") as handle:
                handle.write(a3m_lines)
            # filter
            os.system(f"hhfilter -i {path}/in_{msa_name}_{n}.a3m -id {id} -cov {cov} -o {path}/out_{msa_name}_{n}.a3m")
            with open(f"{path}/out_{msa_name}_{n}.a3m", "r") as handle:
                for line in handle:
                    if not line.startswith(">"):
                        xs = ["-" * l for l in Ls]
                        xs[n] = line.rstrip()
                        # expand homooligomeric sequences
                        xs = ["/".join([x] * num) for x, num in zip(xs, u_nums)]
                        sub_msa[-1].append("/".join(xs))
                        sub_msa_num += 1
            sub_idx.append(list(range(len(sub_msa[-1]))))

        while len(msa) < max_msa and sub_msa_num > 0:
            for n in range(len(sub_idx)):
                if len(sub_idx[n]) > 0:
                    msa.append(sub_msa[n][sub_idx[n].pop(0)])
                    sub_msa_num -= 1
                if len(msa) == max_msa:
                    break

    with open(f"{jobname}/{msa_name}.a3m", "w") as handle:
        for n, sequence in enumerate(msa):
            handle.write(f">n{n}\n{sequence}\n")


def run_mmseqs2(
    x,
    prefix,
    use_env=True,
    use_filter=True,
    use_templates=False,
    filter=None,
    use_pairing=False,
    use_pairwise=False,
    use_taxonomy=False,
    host_url="https://api.colabfold.com",
) -> Tuple[List[str], List[str]]:
    submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"

    def submit(seqs, mode, N=101):
        n, query = N, ""
        for seq in seqs:
            query += f">{n}\n{seq}\n"
            n += 1

        while True:
            error_count = 0
            try:
                # https://requests.readthedocs.io/en/latest/user/advanced/#advanced
                # "good practice to set connect timeouts to slightly larger than a multiple of 3"
                res = requests.post(
                    f"{host_url}/{submission_endpoint}",
                    data={"q": query, "mode": mode},
                    timeout=6.02,
                )
            except requests.exceptions.Timeout:
                logger.warning("Timeout while submitting to MSA server. Retrying...")
                continue
            except Exception as e:
                error_count += 1
                logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break

        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def status(ID):
        while True:
            error_count = 0
            try:
                res = requests.get(f"{host_url}/ticket/{ID}", timeout=6.02)
            except requests.exceptions.Timeout:
                logger.warning("Timeout while fetching status from MSA server. Retrying...")
                continue
            except Exception as e:
                error_count += 1
                logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break
        try:
            out = res.json()
        except ValueError:
            logger.error(f"Server didn't reply with json: {res.text}")
            out = {"status": "ERROR"}
        return out

    def download(ID, path):
        error_count = 0
        while True:
            try:
                res = requests.get(f"{host_url}/result/download/{ID}", timeout=6.02)
            except requests.exceptions.Timeout:
                logger.warning("Timeout while fetching result from MSA server. Retrying...")
                continue
            except Exception as e:
                error_count += 1
                logger.warning(f"Error while fetching result from MSA server. Retrying... ({error_count}/5)")
                logger.warning(f"Error: {e}")
                time.sleep(5)
                if error_count > 5:
                    raise
                continue
            break
        with open(path, "wb") as out:
            out.write(res.content)

    # process input x
    seqs = [x] if isinstance(x, str) else x

    # compatibility to old option
    if filter is not None:
        use_filter = filter

    # setup mode
    if use_filter:
        mode = "env" if use_env else "all"
    else:
        mode = "env-nofilter" if use_env else "nofilter"

    if use_taxonomy:
        mode += "-taxonomy"

    if use_pairing:
        mode = ""
        use_templates = False
        use_env = False

    if use_pairwise:
        mode = "pairwise"
        use_templates = False
        use_env = False

    # define path
    path = f"{prefix}_{mode}"
    if not os.path.isdir(path):
        os.mkdir(path)

    # call mmseqs2 api
    tar_gz_file = f"{path}/out.tar.gz"
    N, REDO = 101, True

    # deduplicate and keep track of order
    seqs_unique = []
    # TODO this might be slow for large sets
    [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
    Ms = [N + seqs_unique.index(seq) for seq in seqs]
    # lets do it!
    if not os.path.isfile(tar_gz_file):
        TIME_ESTIMATE = 150 * len(seqs_unique)
        with tqdm(total=TIME_ESTIMATE, bar_format=TQDM_BAR_FORMAT) as pbar:
            while REDO:
                pbar.set_description("SUBMIT")

                # Resubmit job until it goes through
                out = submit(seqs_unique, mode, N)
                while out["status"] in ["UNKNOWN", "RATELIMIT"]:
                    sleep_time = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
                    # resubmit
                    time.sleep(sleep_time)
                    out = submit(seqs_unique, mode, N)

                if out["status"] == "ERROR":
                    raise Exception(f"MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.")

                if out["status"] == "MAINTENANCE":
                    raise Exception(f"MMseqs2 API is undergoing maintenance. Please try again in a few minutes.")

                # wait for job to finish
                ID, TIME = out["id"], 0
                pbar.set_description(out["status"])
                while out["status"] in ["UNKNOWN", "RUNNING", "PENDING"]:
                    t = 5 + random.randint(0, 5)
                    logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
                    time.sleep(t)
                    out = status(ID)
                    pbar.set_description(out["status"])
                    if out["status"] == "RUNNING":
                        TIME += t
                        pbar.update(n=t)
                    # if TIME > 900 and out["status"] != "COMPLETE":
                    #  # something failed on the server side, need to resubmit
                    #  N += 1
                    #  break

                if out["status"] == "COMPLETE":
                    if TIME < TIME_ESTIMATE:
                        pbar.update(n=(TIME_ESTIMATE - TIME))
                    REDO = False

                if out["status"] == "ERROR":
                    REDO = False
                    raise Exception(f"MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.")

            # Download results
            download(ID, tar_gz_file)

    # prep list of a3m files
    if use_pairing:
        a3m_files = [f"{path}/pair.a3m"]
    else:
        a3m_files = [f"{path}/uniref.a3m"]
        if use_env:
            a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

    # extract a3m files
    if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
        with tarfile.open(tar_gz_file) as tar_gz:
            tar_gz.extractall(path)

    # templates
    if use_templates:
        templates = {}
        # print("seq\tpdb\tcid\tevalue")
        for line in open(f"{path}/pdb70.m8", "r"):
            p = line.rstrip().split()
            M, pdb, qid, e_value = p[0], p[1], p[2], p[10]
            M = int(M)
            if M not in templates:
                templates[M] = []
            templates[M].append(pdb)
            # if len(templates[M]) <= 20:
            #  print(f"{int(M)-N}\t{pdb}\t{qid}\t{e_value}")

        template_paths = {}
        for k, TMPL in templates.items():
            TMPL_PATH = f"{prefix}_{mode}/templates_{k}"
            if not os.path.isdir(TMPL_PATH):
                os.mkdir(TMPL_PATH)
                TMPL_LINE = ",".join(TMPL[:20])
                response = None
                while True:
                    error_count = 0
                    try:
                        # https://requests.readthedocs.io/en/latest/user/advanced/#advanced
                        # "good practice to set connect timeouts to slightly larger than a multiple of 3"
                        response = requests.get(
                            f"{host_url}/template/{TMPL_LINE}",
                            stream=True,
                            timeout=6.02,
                        )
                    except requests.exceptions.Timeout:
                        logger.warning("Timeout while submitting to template server. Retrying...")
                        continue
                    except Exception as e:
                        error_count += 1
                        logger.warning(f"Error while fetching result from template server. Retrying... ({error_count}/5)")
                        logger.warning(f"Error: {e}")
                        time.sleep(5)
                        if error_count > 5:
                            raise
                        continue
                    break
                with tarfile.open(fileobj=response.raw, mode="r|gz") as tar:
                    tar.extractall(path=TMPL_PATH)
                os.symlink("pdb70_a3m.ffindex", f"{TMPL_PATH}/pdb70_cs219.ffindex")
                with open(f"{TMPL_PATH}/pdb70_cs219.ffdata", "w") as f:
                    f.write("")
            template_paths[k] = TMPL_PATH

    # gather a3m lines
    a3m_lines = {}
    for a3m_file in a3m_files:
        update_M, M = True, None
        for line in open(a3m_file, "r"):
            if len(line) > 0:
                if "\x00" in line:
                    line = line.replace("\x00", "")
                    update_M = True
                if line.startswith(">") and update_M:
                    M = int(line[1:].rstrip())
                    update_M = False
                    if M not in a3m_lines:
                        a3m_lines[M] = []
                a3m_lines[M].append(line)

    # return results

    a3m_lines = ["".join(a3m_lines[n]) for n in Ms]

    if use_templates:
        template_paths_ = []
        for n in Ms:
            if n not in template_paths:
                template_paths_.append(None)
                # print(f"{n-N}\tno_templates_found")
            else:
                template_paths_.append(template_paths[n])
        template_paths = template_paths_

    return (a3m_lines, template_paths) if use_templates else a3m_lines
