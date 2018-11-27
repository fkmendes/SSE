import numpy as np
import csv
import argparse
import matplotlib.pyplot as plt
import os

def read(csv_name):
    with open(csv_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                labels = row
            if line_count == 1:
                post = [float(r) for r in row]
            line_count += 1
    return np.array(post), labels

def get_ordering(labels):
    labels = [int(label.replace("nd", "")) for label in labels]
    indices = np.argsort(labels)
    return indices

def prepare_csv(csv_name):  # read and sort the nodes
    post, labels = read(csv_name)
    indices = get_ordering(labels)
    post = post[indices]
    return post

def construct_post_comparison(csv_1_name, csv_2_name, output_file_name):
    csv_1 = prepare_csv(csv_1_name)
    csv_2 = prepare_csv(csv_2_name)

    #for pair in enumerate(zip(csv_1, csv_2)):
    #    print(pair)

    plt.xlim(0, 1)
    plt.ylim(0, 1)

    fig, ax = plt.subplots()
    ax.scatter(csv_1, csv_2, s=25, cmap=plt.cm.coolwarm, zorder=10)

    lims = [
        np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
        np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
    ]

    # now plot both limits against eachother
    ax.plot(lims, lims, 'k-', alpha=0.75, zorder=0)
    ax.set_aspect('equal')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    fig.savefig("plots/" + output_file_name + '.png', dpi=300)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Generate asr validation", description="Generates ancestral state reconstruction figures that validate method.")
    parser.add_argument("-c1", "--csv-1-name", action="store", dest="csv1", default=None, type=str, help="CSV file containing post data.")
    parser.add_argument("-c2", "--csv-2-name", action="store", dest="csv2", default=None, type=str, help="CSV file containing post data.")
    parser.add_argument("-o", "--output-name", action="store", dest="output", default="asr", type=str, help="Name of output file.")
    args = parser.parse_args()

    if not os.path.exists("plots"):
        os.makedirs("plots")

    construct_post_comparison(args.csv1, args.csv2, args.output)
