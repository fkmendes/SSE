import pandas
import matplotlib.pyplot as plt
import argparse
import numpy as np

def plot(name, use_log):
    df = pandas.read_csv(name, index_col=0, header=None)
    n_samples, n_nodes = df.shape
    for i in range(1, n_nodes + 1):
        fig, ax = plt.subplots()
        plt.title(name)
        if use_log:
            ax.scatter(np.log(df.index), df[i])
        else:
            ax.scatter(df.index, df[i])
        fig.savefig("plots/" + name + str(i) + '.png', dpi=300)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Generate asr validation", description="Generates ancestral state reconstruction figures that validate method.")
    parser.add_argument("-n", "--csv-name", action="store", dest="csv", default=None, type=str, help="CSV file containing post data.")
    parser.add_argument("-l", "--log-scale", action="store", dest="log", default=False, type=bool, help=".")
    args = parser.parse_args()

    plot(args.csv, args.log)
