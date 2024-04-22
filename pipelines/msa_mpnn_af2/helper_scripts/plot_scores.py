import matplotlib.pyplot as plt
import seaborn as sns


def create_and_save_plot(df, output):
    colors = ["orange", "lightsteelblue"]
    sns.set_palette(sns.color_palette(colors))

    fig, axes = plt.subplots(1, 3, figsize=(10, 3))
    plt.subplots_adjust(wspace=0.4)

    for ax, metric in zip(axes.flat, ["rmsd", "plddt", "pae"]):
        sns.stripplot(data=df, x="temp", y=metric, ax=ax, dodge=True)

    sns.despine()
    fig.patch.set_facecolor("white")

    # Save the plot
    if output:
        plt.savefig(output, bbox_inches="tight", pad_inches=0.1, dpi=300)
    # Display the plot (optional)
    plt.show()
