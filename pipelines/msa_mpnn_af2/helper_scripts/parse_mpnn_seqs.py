import pandas as pd


def extract_scores(fasta_file, rank=None):
    with open(fasta_file) as f:
        data = f.readlines()

    # Removing the first 2 lines with original sequence
    data = data[2:]
    # Looping through each line to extract the required information
    results = []
    for i in range(0, len(data), 2):
        T = float(data[i].split("T=")[1].split(",")[0])
        sample = int(data[i].split("sample=")[1].split(",")[0])
        score = float(data[i].split("score=")[1].split(",")[0])
        global_score = float(data[i].split("global_score=")[1].split(",")[0])
        seq_recovery = float(data[i].split("seq_recovery=")[1].split("\n")[0])
        sequence = data[i + 1].strip()

        # Adding the extracted values to the list
        results.append([T, sample, score, global_score, seq_recovery, sequence])

    # Converting the list of values to a pandas DataFrame
    columns = ["temp", "sample", "score", "global_score", "sequence_recovery", "seq"]
    df = pd.DataFrame(results, columns=columns)

    # Ranking the DataFrame based on the desired score
    if rank:
        df.sort_values(by=rank, ascending=False, inplace=True)
    df.reset_index(drop=True, inplace=True)
    # df.to_csv(f"{fasta_file.split('.')[0]}_scores.csv", index=False)
    return df


def calculate_seq_identity(seq1, seq2, excluded, exclude_AA="X"):
    count_all = 0
    count_same = 0
    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if i not in excluded and not a == exclude_AA:
            count_all += 1
            if a == b:
                count_same += 1

    if count_all == 0:
        return 0.0

    # (sum(a == b for a, b in zip(seq1, seq2)) / len(seq1))
    return count_same / count_all
